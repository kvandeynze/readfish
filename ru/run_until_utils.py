import logging
import logging.handlers
import os
import toml
import threading
import time
from watchdog.events import FileSystemEventHandler
from ru.mapper import MappingServer as Map
from ru.utils import nice_join, print_args, send_message, Severity, get_device


def file_dict_of_folder_simple(path, args, logging, fastqdict):
    logger = logging.getLogger("ExistingFileProc")

    file_list_dict = dict()

    counter = 0

    if os.path.isdir(path):

        logger.info("caching existing fastq files in: %s" % (path))

        for path, dirs, files in os.walk(path):

            for f in files:

                counter += 1

                if f.endswith(".fastq") or f.endswith(".fastq.gz"):

                    logger.debug("Processing File {}\r".format(f))
                    filepath = os.path.join(path, f)
                    file_list_dict[filepath] = os.stat(filepath).st_mtime

    logger.info("processed %s files" % (counter))



    logger.info("found %d existing fastq files to process first." % (len(file_list_dict)))

    return file_list_dict

def write_new_toml(args,targets):
    for k in args.toml["conditions"].keys():
        curcond = args.toml["conditions"].get(k)
        if isinstance(curcond,dict):

            #newtargets = targets
            #newtargets.extend(curcond["targets"])

            #newtargets = list(dict.fromkeys(newtargets))
            #curcond["targets"]=list(set(newtargets))
            curcond["targets"]=targets

    with open("{}_live".format(args.tomlfile), "w") as f:
        toml.dump(args.toml,f)



###Function modified from https://raw.githubusercontent.com/lh3/readfq/master/readfq.py


def readfq(fp):  # this is a generator function
    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for l in fp:  # search for the start of the next record
                if l[0] in ">@":  # fasta/q header line
                    last = l[:-1]  # save this line
                    break
        if not last:
            break
        desc, name, seqs, last = last[1:], last[1:].partition(" ")[0], [], None
        for l in fp:  # read the sequence
            if l[0] in "@+>":
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != "+":  # this is a fasta record
            yield desc, name, "".join(seqs), None  # yield a fasta record
            if not last:
                break
        else:  # this is a fastq record
            seq, leng, seqs = "".join(seqs), 0, []
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield desc, name, seq, "".join(seqs)  # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield desc, name, seq, None  # yield a fasta record instead
                break


def fastq_results(fastq):
    if fastq.endswith(".gz"):

        with gzip.open(fastq, "rt") as fp:
            try:
                for desc, name, seq, qual in readfq(fp):
                    yield desc, name, seq, qual

            except Exception as e:
                print(e)
    else:
        with open(fastq, "r") as fp:
            try:
                # now = time.time()
                for desc, name, seq, qual in readfq(fp):
                    yield desc, name, seq, qual

            except Exception as e:
                print(e)


def parse_fastq_file(fastqfilelist,args,logging,mapper):
    logger = logging.getLogger("ParseFastq")
    # Add the reference to the mapper
    #ToDo: This needs to be some kind of real reference name.
    mapper.add_reference("test",args.toml['conditions']['reference'])

    for file in fastqfilelist:
        for desc,name, seq,qual in fastq_results(file):
            sequence_list=({"sequence":seq,"read_id":name})
            mapper.map_sequence("test",sequence_list)

class FastqHandler(FileSystemEventHandler):

    def __init__(self, args,logging,rpc_connection):
        self.args = args
        #self.messageport = messageport
        self.connection = rpc_connection
        self.logger = logging.getLogger("FastqHandler")
        self.running = True
        self.fastqdict = dict()
        self.mapper=Map()
        self.mapper.set_cov_target(args.depth)
        self.creates = file_dict_of_folder_simple(self.args.watch, self.args, logging,
                                                  self.fastqdict)
        self.t = threading.Thread(target=self.processfiles)

        try:
            self.t.start()
        except KeyboardInterrupt:
            self.t.stop()
            raise

    def on_created(self, event):
        """Watchdog counts a new file in a folder it is watching as a new file"""
        """This will add a file which is added to the watchfolder to the creates and the info file."""
        # if (event.src_path.endswith(".fastq") or event.src_path.endswith(".fastq.gz")):
        #     self.creates[event.src_path] = time.time()


        # time.sleep(5)
        if (event.src_path.endswith(".fastq") or event.src_path.endswith(".fastq.gz") or event.src_path.endswith(
                ".fq") or event.src_path.endswith(".fq.gz")):
            self.logger.info("Processing file {}".format(event.src_path))
            self.creates[event.src_path] = time.time()

    def on_modified(self, event):
        if (event.src_path.endswith(".fastq") or event.src_path.endswith(".fastq.gz") or event.src_path.endswith(
                ".fq") or event.src_path.endswith(".fq.gz")):
            self.logger.info("Processing file {}".format(event.src_path))
            self.logger.debug("Modified file {}".format(event.src_path))
            self.creates[event.src_path] = time.time()

    def on_moved(self, event):
        if any((event.dest_path.endswith(".fastq"), event.dest_path.endswith(".fastq,gz"),
                event.dest_path.endswith(".gq"), event.dest_path.endswith(".fq.gz"))):
            self.logger.info("Processing file {}".format(event.dest_path))
            self.logger.debug("Modified file {}".format(event.dest_path))
            self.creates[event.dest_path] = time.time()