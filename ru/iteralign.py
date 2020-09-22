"""
Iter-align.
Design spec.
    1. Grab fastq files from a location.
    2. Align the files ONCE to a reference.
    3. Calculate cumulative coverage information.
    4. Write a list of genomes that are covered at a particular threshold.
    5. Rinse and repeat
"""
import logging
import logging.handlers
import os
import sys
import time
import toml
from watchdog.observers.polling import PollingObserver as Observer

from ru.utils import nice_join, print_args, send_message, Severity, get_device
from minknow_api.acquisition_pb2 import MinknowStatus
from ru.run_until_utils import parse_fastq_file,write_new_toml
from ru.run_until_utils import FastqHandler
from ru.ru_gen import _cli as BASE
from ru.ru_gen import run as dnrun
from argparse import Namespace


DEFAULT_SERVER_HOST = "127.0.0.1"
DEFAULT_LOG_FORMAT = "%(asctime)s - %(name)-20s - %(message)s"
LOG_LEVELS = ("debug", "info", "warning", "error", "critical")
DEFAULT_COVERAGE_DEPTH = 30
DEFAULT_PERCENTAGE_COVERED = 0.99
DEFAULT_CORES = 2


_help = "ReadFish and Run Until, using minimap2"
_cli = BASE + (

    #(
    #    "--host",
    #    dict(
    #        metavar="HOST",
    #        help="MinKNOW server host",
    #        default=DEFAULT_SERVER_HOST,
    #    ),
    #),
    #(
    #    "--device",
    #    dict(
    #        metavar="DEVICE",
    #        type=str,
    #        help="Name of the sequencing position e.g. MS29042 or GA10000 etc.",
    #        required=True,
    #    ),
    #),

    (
        "--watch",
        dict(
            metavar="FOLDER",
            help="Top Level Folder containing fastq reads.",
            default=None,
        ),
    ),
    (
        "--percent",
        dict(
            metavar="PERCENT",
            help="Default percent of target covered at given depth (default {})".format(DEFAULT_PERCENTAGE_COVERED),
            default=DEFAULT_PERCENTAGE_COVERED,
            type=float,
        ),
    ),
    (
        "--depth",
        dict(
            metavar="DEPTH",
            help="Desired coverage depth (default {})".format(DEFAULT_COVERAGE_DEPTH),
            default=DEFAULT_COVERAGE_DEPTH,
            type=int,
        ),
    ),
    (
        "--threads",
        dict(
            metavar="THREADS",
            help="Set the number of default threads to use for threaded tasks (default {})".format(DEFAULT_CORES),
            default=DEFAULT_CORES,
            type=int,
        ),
    ),
    #(
    #    "--log-level",
    #    dict(
    #        metavar="LOG-LEVEL",
    #        action="store",
    #        default="info",
    #        choices=LOG_LEVELS,
    #        help="One of: {}".format(nice_join(LOG_LEVELS)),
    #    ),
    #),
    #(
    #    "--log-format",
    #    dict(
    #        metavar="LOG-FORMAT",
    #        action="store",
    #        default=DEFAULT_LOG_FORMAT,
    #        help="A standard Python logging format string (default: {!r})".format(
    #            DEFAULT_LOG_FORMAT.replace("%", "%%")
    #        ),
    #    ),
    #),
    #(
    #    "--log-file",
    #    dict(
    #        metavar="LOG-FILE",
    #        action="store",
    #        default=None,
    #        help="A filename to write logs to, or None to write to the standard stream (default: None)",
    #    ),
    #),
    #(
    #    "--toml",
    #    dict(
    #        metavar="TOML",
    #        required=True,
    #        help="The magic TOML file that will save your life?",
    #        #type=toml.load,
    #    ),
    #),
)


class FastQMonitor(FastqHandler):
    def __init__(self,args,logging,rpc_connection):
        super().__init__(
            args=args,
            logging=logging,
            #messageport=messageport,
            rpc_connection=rpc_connection,
        )

    def processfiles(self):
        self.logger.info("Process Files Initiated")
        self.counter = 0
        self.targets = []

        while self.running:
            currenttime = time.time()
            fastqfilelist=list()
            for fastqfile, createtime in sorted(self.creates.items(), key=lambda x: x[1]):
                delaytime = 0
                # file created 5 sec ago, so should be complete. For simulations we make the time longer.
                if (int(createtime) + delaytime < time.time()):
                    self.logger.info(fastqfile)
                    del self.creates[fastqfile]
                    self.counter +=1
                    fastqfilelist.append(fastqfile)
            parse_fastq_file(fastqfilelist,self.args,logging,self.mapper)
            #print (self.mapper.report_coverage())
            #This prints those targets with a coverage greater than the threshold set in the arguments
            targets = self.mapper.target_coverage().keys()
            if len(targets) > len(self.targets):
                updated_targets = set(targets) - set(self.targets)
                update_message = "Updating targets with {}".format(nice_join(updated_targets, conjunction="and"))
                self.logger.info(update_message)
                if not self.args.simulation:
                    send_message(self.connection, update_message, Severity.WARN)
                write_new_toml(self.args,targets)
                self.targets = []
                self.targets = targets

            if len(self.targets)>0 and self.mapper.check_complete():
                self.logger.info("Every target is covered at at least {}x".format(self.args.depth))
                if not self.args.simulation:
                    self.connection.protocol.stop_protocol()
                    send_message(
                        self.connection,
                        "Iter Align has stopped the run as all targets should be covered by at least {}x".format(
                            self.args.depth
                        ),
                        Severity.WARN,
                    )

            if currenttime+5 > time.time():
                time.sleep(5)



def main():
    sys.exit(
        "This entry point is deprecated, please use 'readfish align' instead"
    )


def run(parser, args):
    args_copy = Namespace(**vars(args))
    args.tomlfile = args.toml
    args.toml = toml.load(args.toml)
    print (args)

    # TODO: Move logging config to separate configuration file
    # set up logging to file
    logging.basicConfig(level=logging.DEBUG,
                        format='%(levelname)s::%(asctime)s::%(name)s::%(message)s',
                        filename=args.log_file,
                        filemode='w')

    # define a Handler that writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)

    # set a format which is simpler for console use
    formatter = logging.Formatter('%(name)-15s: %(levelname)-8s %(message)s')
    console.setFormatter(formatter)

    # add the handler to the root logger
    logging.getLogger('').addHandler(console)

    # Start by logging sys.argv and the parameters used
    logger = logging.getLogger("Manager")
    logger.info(" ".join(sys.argv))

    logger.info("Initialising iterAlign.")

    logger.info("Setting up FastQ monitoring.")

    #### Check if a run is active - if not, wait.

    args.simulation = True
    connection = None
    if args.watch is None:
        args.simulation = False
        logger.info("Creating rpc connection for device {}.".format(args.device))
        try:
            connection = get_device(args.device).connect()
        except ValueError as e:
            print(e)
            sys.exit(1)

        send_message(connection, "Iteralign Connected to MinKNOW", Severity.WARN)

        logger.info("Loaded RPC")
        while connection.acquisition.current_status().status != MinknowStatus.PROCESSING:
            time.sleep(1)
        #### Check if we know where data is being written to , if not... wait
        args.watch = connection.acquisition.get_acquisition_info().config_summary.reads_directory



    event_handler = FastQMonitor(args,logging,connection)
    # This block handles the fastq
    observer = Observer()
    print (args.watch)
    observer.schedule(event_handler, path=args.watch, recursive=True)
    observer.daemon = True




    try:

        observer.start()
        logger.info("FastQ Monitoring Running.")

        dnrun(parser, args_copy)

        while 1:
            time.sleep(1)

    except KeyboardInterrupt:

        logger.info("Exiting - Will take a few seconds to clean up!")

        observer.stop()
        observer.join()

        os._exit(0)


