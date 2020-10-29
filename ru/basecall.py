"""basecall.py

Extension of pyguppy Caller that maintains a connection to the basecaller

"""
import os
import logging
import time
from collections import namedtuple
from multiprocessing import Process, Queue
import queue


import mappy as mp
import numpy as np

try:
    from pyguppy_client_lib.pyclient import PyGuppyClient
    from pyguppy_client_lib.helper_functions import package_read

    pyguppy_available = True
except ImportError:
    pyguppy_available = False


try:
    import deepnano2

    deepnano2_available = True
except ImportError:
    deepnano2_available = False


__all__ = ["GuppyCaller", "DeepNanoCaller"]

logger = logging.getLogger("RU_basecaller")
CALIBRATION = namedtuple("calibration", "scaling offset")


class DefaultDAQValues:
    """Provides default calibration values

    Mimics the read_until_api calibration dict value from
    https://github.com/nanoporetech/read_until_api/blob/2319bbe80889a17c4b38dc9cdb45b59558232a7e/read_until/base.py#L34
    all keys return scaling=1.0 and offset=0.0
    """

    calibration = CALIBRATION(1.0, 0.0)

    def __getitem__(self, _):
        return self.calibration


def _call(queue_in, queue_out, net_type, path, beam, cut):
    caller = deepnano2.Caller(net_type, path, beam, cut)
    while True:
        try:
            read = queue_in.get()
        except queue.Empty:
            continue
        read["seq"], read["qual"] = caller.call_raw_signal(rescale_signal(read["signal"]))
        queue_out.put(read)


def med_mad(x, factor=1.4826):
    """
    Calculate signal median and median absolute deviation
    """
    med = np.median(x)
    mad = np.median(np.absolute(x - med)) * factor
    return med, mad


def rescale_signal(signal):
    med, mad = med_mad(signal)
    signal -= med
    signal /= mad
    return signal


class Caller:
    def basecall_minknow(self, reads, signal_dtype, decided_reads, **kwargs):
        raise NotImplementedError


class DeepNanoCaller(Caller):
    def __init__(
        self, network_type="56", beam_size=5, beam_cut_threshold=0.01, threads=2
    ):
        if not deepnano2_available:
            # Could use sys.exit here instead
            raise ImportError("DeepNano2 not available")
        self.net_type = network_type
        self.beam = beam_size
        self.cut = beam_cut_threshold
        # TODO: do this better
        self.net_path = os.path.join(
            deepnano2.__path__[0], "weights", "rnn{}.txt".format(self.net_type)
        )
        self.threads = threads
        # TODO: add max queue size for _in and out
        self._in, self._out = Queue(), Queue()
        args = (self._in, self._out, self.net_type, self.net_path, self.beam, self.cut)
        self.workers = [Process(target=_call, args=args) for _ in range(self.threads)]

    def __repr__(self):
        return "{}()".format(self.__class__.__name__)

    def start(self):
        for worker in self.workers:
            worker.start()

    def terminate(self):
        for worker in self.workers:
            worker.terminate()

    def pass_data(self, reads):
        """Put reads into the queue for calling

        Parameters
        ----------
        reads : Iterable[Dict]

        Returns
        -------
        count : int
            The number of reads put in the queue
        """
        count = 0
        for r in reads:
            self._in.put(r)
            count += 1
        return count

    def get_completed(self):
        # Get a list of dict back with seq and qual added
        items = []
        while True:
            try:
                items.append(self._out.get_nowait())
            except queue.Empty:
                break
        return items

    @staticmethod
    def pack_read(reads, signal_dtype, decided_reads):
        """
        Parameters
        ----------
        reads
        signal_dtype
        Returns
        -------
        """
        for channel, read in reads:
            if decided_reads.get(channel, "") == read.id:
                continue
            yield {
                "signal": np.frombuffer(read.raw_data, dtype=signal_dtype),
                "read_id": read.id,
                "channel": channel,
                "read_number": read.number,
            }

    def basecall_minknow(self, reads, signal_dtype, decided_reads, **kwargs):
        """DeepNano basecaller wrapper for MinKNOW RPC reads

        Parameters
        ----------
        **kwargs
        reads : iterable[Tuple[int, rpc.Read]]
            List or generator of tuples containing (channel, MinKNOW.rpc.Read)
        signal_dtype
            Numpy dtype of the raw data
        prev_signal : DefaultDict[int: collections.deque[Tuple[str, np.ndarray]]]
            Dictionary of previous signal fragment from a channel
        decided_reads : Dict[int: str]
            Dictionary of channels with the last read id a decision was made for
        Yields
        ------
        read_info : tuple
            Tuple of read info (channel, read_number)
        read_id : str
        sequence : str
        sequence_length : int
        quality : str
        """
        done = 0
        sent = self.pass_data(self.pack_read(reads, signal_dtype, decided_reads))

        while done < sent:
            res = self.get_completed()

            if not res:
                continue
            for read in res:
                yield (read["channel"], read["read_number"]), read["read_id"], read[
                    "seq"
                ], len(read["seq"]), read["qual"]
                done += 1

    def disconnect(self):
        self.terminate()


class GuppyCaller(PyGuppyClient):
    # TODO: Use Caller class as well
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.connect()

    def basecall_minknow(self, reads, signal_dtype, decided_reads, **kwargs):
        """Guppy basecaller wrapper for MinKNOW RPC reads

        Parameters
        ----------
        reads : iterable[Tuple[int, rpc.Read]]
            List or generator of tuples containing (channel, MinKNOW.rpc.Read)
        signal_dtype
            Numpy dtype of the raw data
        decided_reads : Dict[int: str]
            Dictionary of channels with the last read id a decision was made for

        Yields
        ------
        read_info : tuple
            Tuple of read info (channel, read_number)
        read_id : str
        sequence : str
        sequence_length : int
        quality : str
        """
        hold = {}
        # FixMe: This is resolved in later versions of guppy.
        skipped = {}
        done = 0
        read_counter = 0

        daq_values = kwargs.get("daq_values", DefaultDAQValues())

        for channel, read in reads:
            hold[read.id] = (channel, read.number)
            t0 = time.time()
            success = self.pass_read(
                package_read(
                    read_id=read.id,
                    raw_data=np.frombuffer(read.raw_data, signal_dtype),
                    daq_offset=daq_values[channel].offset,
                    daq_scaling=daq_values[channel].scaling,
                )
            )
            if not success:
                logging.warning("Skipped a read: {}".format(read.id))
                # FixMe: This is resolved in later versions of guppy.
                skipped[read.id] = hold.pop(read.id)
                continue
            else:
                read_counter += 1

            sleep_time = self.throttle - t0
            if sleep_time > 0:
                time.sleep(sleep_time)

        while done < read_counter:
            results = self.get_completed_reads()

            if not results:
                time.sleep(self.throttle)
                continue

            for r in results:
                r_id = r["metadata"]["read_id"]
                try:
                    i = hold.pop(r_id)
                except KeyError:
                    # FixMe: This is resolved in later versions of guppy.
                    i = skipped.pop(r_id)
                    read_counter += 1
                yield (
                    i,
                    r_id,
                    r["datasets"]["sequence"],
                    r["metadata"]["sequence_length"],
                    r["datasets"]["qstring"],
                )
                done += 1


class Mapper:
    def __init__(self, index):
        self.index = index
        if self.index:
            self.mapper = mp.Aligner(self.index, preset="map-ont")
            self.initialised = True
        else:
            self.mapper = None
            self.initialised = False

    def map_read(self, seq):
        return self.mapper.map(seq)

    def map_reads(self, calls):
        for read_id, seq in calls:
            yield read_id, list(self.mapper.map(seq))

    def map_reads_2(self, calls):
        """Align reads against a reference

        Parameters
        ----------
        calls : iterable [tuple,  str, str, int, str]
            An iterable of called reads from PerpetualCaller.basecall_minknow

        Yields
        ------
        read_info : tuple
            Tuple of read info (channel, read_number)
        read_id : str
        sequence : str
        sequence_length : int
        mapping_results : list
        """
        for read_info, read_id, seq, seq_len, quality in calls:
            yield read_info, read_id, seq_len, list(self.mapper.map(seq))
