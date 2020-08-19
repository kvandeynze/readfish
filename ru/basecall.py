"""basecall.py

Extension of pyguppy Caller that maintains a connection to the basecaller

"""
import logging
import os

import mappy as mp
import numpy as np

from pyguppyclient.client import GuppyBasecallerClient
from pyguppyclient.decode import ReadData as GuppyRead
from ru.dnb_mp import DeepNanoClient


__all__ = ["GuppyCaller"]

logger = logging.getLogger("RU_basecaller")


def _create_deepnano_read(reads, signal_dtype,signalchunk=4000):
    """

    Parameters
    ----------
    reads
    signal_dtype

    Returns
    -------

    """
    for channel,read in reads:
        yield {
            "signal":np.frombuffer(read.raw_data, dtype=signal_dtype)[-signalchunk:],
            "read_id":read.id,
            "channel":channel,
            "read_number":read.number,
        }

def _create_guppy_read(reads, signal_dtype):
    """Convert a read from MinKNOW RPC to GuppyRead

    Parameters
    ----------
    reads : List[Tuple[int, minknow.rpc.read]]
        List of Tuple, containing (channel, read)
    signal_dtype : np.dtype
        A dtype that can be used by numpy to convert the raw data
    previous_signal : dict
        Dict containing previous signal segments

    Yields
    ------
    channel : int
    read_number : int
    GuppyRead
    """
    for channel, read in reads:
        read_obj = GuppyRead(np.frombuffer(read.raw_data, dtype=signal_dtype), read.id, 0, 1)
        yield channel, read.number, read_obj

def _concat_signal(reads, signal_dtype, previous_signal):
    for channel, read in reads:
        old_read_id, old_signal = previous_signal.get(channel, (("", np.empty(0, dtype=signal_dtype)),))[0]

        if old_read_id == read.id:
            signal = np.concatenate((old_signal, np.frombuffer(read.raw_data, dtype=signal_dtype)))
        else:
            signal = np.frombuffer(read.raw_data, dtype=signal_dtype)

        yield read.id, channel, read.number, signal

def _process_signal(reads, signal_dtype):
    for channel, read in reads:
        signal = np.frombuffer(read.raw_data, dtype=signal_dtype)
        yield read.id, channel, read.number, signal[-4000:]

def _trim_blank(sig, window=300):
    N = len(sig)
    variances = [np.var(sig[i:i+window]) for i in range(N//2, N-window, window)]
    mean_var = np.mean(variances)
    trim_idx = 20
    while window > 5:
        while np.var(sig[trim_idx: trim_idx + window]) < 0.3*mean_var:
            trim_idx += 1
        window //= 2

    return trim_idx


def _trim(signal, window_size=40, threshold_factor=3.0, min_elements=3):
    med, mad = _med_mad(signal[-(window_size*25):])
    threshold = med + mad * threshold_factor
    num_windows = len(signal) // window_size

    for pos in range(num_windows):

        start = pos * window_size
        end = start + window_size

        window = signal[start:end]

        if len(window[window > threshold]) > min_elements:
            if window[-1] > threshold:
                continue
            return end

    return 0


def _rescale_signal(signal):
    """Rescale signal for DeepNano"""
    signal = signal.astype(np.float32)
    med, mad = _med_mad(signal)
    signal -= med
    signal /= mad
    return signal


def _med_mad(x, factor=1.4826):
    """Calculate signal median and median absolute deviation"""
    med = np.median(x)
    mad = np.median(np.absolute(x - med)) * factor
    return med, mad

class DeepCaller(DeepNanoClient):
    def __init__(
            self,
            network_type="96",
            beam_size=5,
            beam_cut_threshold=0.01,
            threads = 4,
    ):
        super().__init__(
            network_type=network_type,
            beam_size=beam_size,
            beam_cut_threshold=beam_cut_threshold,
            threads=threads
                         )
        self.start()

    def basecall_minknow(self, reads, signal_dtype, decided_reads):
        """DeepNano basecaller wrapper for MinKNOW RPC reads

        Parameters
        ----------
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
        read_counter = 0
        read_list = []
        for read in _create_deepnano_read(reads, signal_dtype):
            if read["read_id"] == decided_reads.get(read["channel"], ""):
                continue
            read_list.append(read)

        self.pass_data(read_list)

        read_counter = len(read_list)

        while done < read_counter:
            res = self.get_completed()

            if not res:
                continue
            for read in res:
                yield (read["channel"],read["read_number"]),read["read_id"],read["seq"],len(read["seq"]),read["qual"]
                done += 1

    def disconnect(self):
        self.terminate()


class CPUCaller():
    import deepnano2

    def __init__(
            self,
            network_type="56",
            beam_size=1,
            beam_cut_threshold=0.01,
    ):
        self.network_type = network_type
        self.beam_size = beam_size
        self.beam_cut_threshold = beam_cut_threshold
        self.caller = self.deepnano2.Caller(
            self.network_type,
            os.path.join(self.deepnano2.__path__[0], "weights", "rnn%s.txt" % self.network_type),
            self.beam_size,
            self.beam_cut_threshold,
        )
        logger.info("CPU Caller Up")

    def basecall_minknow(self, reads, signal_dtype, decided_reads):
        hold = {}
        for read_id, channel, read_number, signal in _process_signal(reads, signal_dtype):
            if read_id == decided_reads.get(channel, ""):
                continue

            start = _trim(signal)
            signal = _rescale_signal(signal)
            #print (len(signal))

            hold[read_id] = (channel, read_number)
            seq = self.caller.call_raw_signal(signal)

            yield hold.pop(read_id), read_id, seq, len(seq), ""

    def disconnect(self):
        """Pass through to make CPU caller compatible with GPU"""
        pass

class GuppyCaller(GuppyBasecallerClient):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.connect()

    def basecall_minknow(self, reads, signal_dtype, decided_reads):
        """Guppy basecaller wrapper for MinKNOW RPC reads

        Parameters
        ----------
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
        read_counter = 0

        hold = {}
        for channel, read_number, read in _create_guppy_read(reads, signal_dtype):
            if read.read_id == decided_reads.get(channel, ""):
                continue

            hold[read.read_id] = (channel, read_number)
            try:
                self.pass_read(read)
            except Exception as e:
                logger.warning("Skipping read: {} due to {}".format(read.read_id, e))
                hold.pop(read.read_id)
                continue
            read_counter += 1

        while done < read_counter:
            res = self._get_called_read()

            if res is None:
                continue

            read, called = res

            yield hold.pop(
                read.read_id
            ), read.read_id, called.seq, called.seqlen, called.qual
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
