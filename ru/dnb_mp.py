from multiprocessing import Process, Queue
import queue
import logging
import os, sys
import time
from uuid import uuid4

import numpy as np
import deepnano2


logging.basicConfig(
    format="%(asctime)s:%(processName)-11s:%(message)s", level=logging.INFO
)
LOG = logging.getLogger()


def _call(_in, out, n, p, bs, bct):
    caller = deepnano2.Caller(n, p, bs, bct,)
    while True:
        try:
            pack = _in.get()
        except queue.Empty:
            continue
        pack["seq"], pack["qual"] = caller.call_raw_signal(rescale_signal(pack["signal"]))
        out.put(pack)


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


class CPUCaller:
    def __init__(
        self, network_type="56", beam_size=5, beam_cut_threshold=0.01, threads=2
    ):
        self.net_type = network_type
        self.beam = beam_size
        self.cut = beam_cut_threshold
        self.net_path = os.path.join(
            deepnano2.__path__[0], "weights", f"rnn{self.net_type}.txt"
        )
        self.threads = threads
        # TODO: add max queue size for _in and out
        self._in, self.out = Queue(), Queue()
        args = (self._in, self.out, self.net_type, self.net_path, self.beam, self.cut)
        self.workers = [Process(target=_call, args=args) for _ in range(self.threads)]

    def __repr__(self):
        return f"{self.__class__.__name__}()"

    def start(self):
        for worker in self.workers:
            worker.start()

    def terminate(self):
        for worker in self.workers:
            worker.terminate()

    def pass_data(self, reads):
        # Pass in a dict with the signal and metadata
        for r in reads:
            self._in.put(r)

    def get_completed(self):
        # Get a list of dict back with seq and qual added
        items = []
        while True:
            try:
                items.append(self.out.get_nowait())
            except queue.Empty:
                break
        return items


def gen(n=1):
    for i in range(n):
        yield {
            "signal": np.random.normal(
                size=4000,
            ).astype("float32"),
            "number": i,
            "read_id": str(uuid4()),
        }


if __name__ == "__main__":
    t = 1
    if len(sys.argv) > 1:
        t = int(sys.argv[1])
    # TODO: test with actual signal data
    c = CPUCaller(threads=t)
    c.start()

    for i in range(4):
        s = 0
        n = 7 ** i
        g = list(gen(n))
        t0 = time.time()
        c.pass_data(g)
        while s < n:
            r = c.get_completed()
            s += len(r)
        LOG.info(f"{time.time() - t0:.8f}s for {n} reads")
    c.terminate()
