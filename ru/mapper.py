import mappy as mp
import time
import threading
from pathlib import Path
import toml

from ru.utils import (
    print_args,
    get_run_info,
    between,
    setup_logger,
    describe_experiment,
)


# a = mp.Aligner("test/MT-human.fa")  # load or build index
# if not a: raise Exception("ERROR: failed to load/build index")
# s = a.seq("MT_human", 100, 200)     # retrieve a subsequence from the index
# print(mp.revcomp(s))                # reverse complement
# for name, seq, qual in mp.fastx_read("test/MT-orang.fa"): # read a fasta/q sequence
#        for hit in a.map(seq): # traverse alignments
#                print("{}\t{}\t{}\t{}".format(hit.ctg, hit.r_st, hit.r_en, hit.cigar_str))


class MappingServer:
    """
    This class exists to provide a stable list of references to map against.
    The class should exist in perpetuity and ensure that we only ever have one instance of a reference loaded in memory.
    """

    def __init__(self):
        self.references = set()  # An index of all references available
        self.coverage = (
            dict()
        )  # A holder which will track coverage of each sequence in a reference
        self.targets = dict()
        self.mappingobjects = dict()
        self.interval = 1800  # Check every 30 minutes to see if a reference has been used - delete it if it hasn't been used.
        self.cov_target = 0
        self.conditions = None
        self.run_info = None

        databasemonitor = threading.Thread(target=self.databasemonitor, args=())
        databasemonitor.daemon = True  # Daemonize thread
        databasemonitor.start()
        # self.references.add("camel")

    def parse_targets(self,tomlfile):
        self.targets_to_skip = dict()
        self.run_info, self.conditions, self.reference, self.caller_kwargs = get_run_info(
            tomlfile, num_channels=512
        )
        self.live_toml = Path("{}_live".format(tomlfile))
        self.tomlfile = tomlfile
        print (self.conditions)
        target_list = []
        for condition in self.conditions:
            print (condition)
            print (condition.coords)
            for strand in condition.coords:
                for chr in condition.coords[strand]:
                    for coord in condition.coords[strand][chr]:
                        print ("{},{},{},{}".format(chr,coord[0],coord[1],strand))
                        target_list.append("{},{},{},{}".format(chr,coord[0],coord[1],strand))

        self.write_new_toml(target_list,self.tomlfile)
        #do we care about strand for tracking coverage do we care about orientation.
        #for

    def write_new_toml(self,targets,tomlfile):
        new_toml = toml.load(tomlfile)
        for k in new_toml["conditions"].keys():
            curcond = new_toml["conditions"].get(k)
            if isinstance(curcond, dict):
                # newtargets = targets
                # newtargets.extend(curcond["targets"])

                # newtargets = list(dict.fromkeys(newtargets))
                # curcond["targets"]=list(set(newtargets))
                curcond["targets"] = targets

        with open(self.live_toml, "w") as f:
            toml.dump(new_toml, f)


    def set_cov_target(self, coverage):
        self.cov_target = coverage

    def get_cov_target(self):
        return self.cov_target

    def check_complete(self):
        return len(self.coverage) == len(self.target_coverage().keys())

    def databasemonitor(self):
        while True:
            print("checking references")
            poplist = list()
            print(self.references)
            for reference in self.mappingobjects:
                print(reference)
                if (
                    self.mappingobjects[reference]["last_used"]
                    < time.time() - self.interval
                ):
                    print("This reference is old.")
                    poplist.append(reference)
                else:
                    print("This reference is OK.")
            for reference in poplist:
                # Delete order important here. Removing from the reference list prevents you trying to map something as it is removed.
                self.delete_reference(reference)
                self.mappingobjects.pop(reference)

            time.sleep(self.interval)

    def add_reference(self, reference, filepath):
        """
        Add a reference to the available set of references.
        :param reference: a string demarking a reference
        :return:
        """
        if reference not in self.references:
            self.references.add(reference)
            self.load_index(reference, str(filepath))

    def add_ref_coverage(self, refname, reflen):
        """
        Adds a specific reference name to the coverage dictionary if it isn't already there.
        This assumes that only one user is connecting to the system and tracking coverage.
        :param refname: name of reference being added to the dictionary
        :return:
        """
        if refname not in self.coverage:
            self.coverage[refname] = dict()
            self.coverage[refname]["bases"] = 0
            self.coverage[refname]["length"] = reflen

    def increment_ref_coverage(self, refname, maplen):
        """
        Updates the coverage dictionary for a specific reference element
        :param refname: the index name for the sequence
        :param reflen: the length of the sequence
        :param maplen: the length of the mapping
        :return:
        """
        self.coverage[refname]["bases"] += maplen
        # self.coverage[refname]["length"]=reflen

    def report_coverage(self):
        coverage_results = []
        for refname in self.coverage:
            if self.coverage[refname]["length"] > 0:
                coverage_results.append(
                    {
                        "refname": refname,
                        "bases": self.coverage[refname]["bases"],
                        "coverage": self.coverage[refname]["bases"]
                        / self.coverage[refname]["length"],
                    }
                )
        if self.conditions:
            print (self.targets)
        return coverage_results

    def target_coverage(self):
        """
        Return targets covered by at least the targetvalue
        Parameters
        ----------
        targetvalue

        Returns
        -------

        """
        coverage_results = dict()
        for refname in self.coverage:
            if self.coverage[refname]["length"] > 0:
                if (
                    self.coverage[refname]["bases"] / self.coverage[refname]["length"]
                    >= self.cov_target
                ):
                    coverage_results[refname] = (
                        self.coverage[refname]["bases"]
                        / self.coverage[refname]["length"]
                    )
        return coverage_results

    def delete_reference(self, reference):
        """
        Remove a reference from the available set of references.
        :param reference: a string demarking a reference
        :return:
        """
        self.references.remove(reference)

    def list_references(self):
        return self.references

    def valid(self, ref_name):
        if ref_name in self.references:
            return True

    def load_index(self, reference, filepath):
        self.mappingobjects[reference] = dict()
        self.mappingobjects[reference]["reference"] = mp.Aligner(
            filepath,
            preset="map-ont",
        )
        self.mappingobjects[reference]["last_used"] = time.time()
        for refname in self.mappingobjects[reference]["reference"].seq_names:
            self.add_ref_coverage(
                refname, len(self.mappingobjects[reference]["reference"].seq(refname))
            )
            # print (len(self.mappingobjects[reference]["reference"].seq(refname)))
            # print (refname)
        # print (self.coverage)
        # print (len(self.coverage))

    def check_tracking(self):
        if self.conditions:
            return True

    def map_sequence(self, reference, sequence, trackcov=True):
        """
        This is a fast mapper that takes a sequence and returns the mapped sequence.
        :param reference: index into the reference dictionary
        :param sequence: list of sequences
        :return: list of map objects
        """
        self.refresh_index(reference)
        results = list()
        for hit in self.mappingobjects[reference]["reference"].map(
            sequence["sequence"]
        ):
            results.append(
                "{}\t{}\t{}".format(sequence["read_id"], len(sequence["sequence"]), hit)
            )
            if trackcov:
                #### How do we handle multiple hits?
                # print ("updating coverage {} {}".format(hit.ctg,hit.mlen))
                self.increment_ref_coverage(hit.ctg, hit.mlen)
                if self.conditions:
                    for item in sequence["desc"].split():
                        if item.startswith('ch='):
                            channel = int(item.split("=")[1])
                            condition = self.run_info[channel]
                            result = False
                            #print (condition,hit.strand,hit.ctg,hit.r_st,hit.r_en,hit.mlen)
                            if hit.strand == 1:
                                result = self.check_target('+',hit.ctg,hit.mlen,hit.r_st,condition)
                                #if result:
                                #    print ("forward match found")
                                pass
                            else:
                                result = self.check_target('-', hit.ctg, hit.mlen, hit.r_en, condition)
                                #if result:
                                #    print ("reverse match found")
                                pass
                            if not result:
                                #print ("no mathc found")
                                pass

        # print(self.report_coverage())
        return results

    def check_target(self,strand,ctg,length,coord,cond):
        """coord_match = any(
            between(int(coord), c)
                for c in self.conditions[cond]
                .coords.get(strand, {})
                .get(ctg, [])
        )"""
        coord_match = [
            c
            for c in self.conditions[cond].coords.get(strand, {}).get(ctg, [])
            if between(int(coord), c)
        ]

        if coord_match:

            left_coord = coord_match[0][0]
            target_length = coord_match[0][1]-coord_match[0][0]+1

            #Going to ignore strand - this assumes the targets overlap. For now we will keep it in.
            #ToDo: properly handle stranded information. For now just ignoring.
            #if strand not in self.targets.keys():
            #    self.targets[strand]=dict()
            if ctg not in self.targets.keys():
                self.targets[ctg]=dict()
            if left_coord not in self.targets[ctg].keys():
                self.targets[ctg][left_coord]=dict()
                self.targets[ctg][left_coord]["cov"]=0
                self.targets[ctg][left_coord]["size"]=target_length
            self.targets[ctg][left_coord]["cov"]+=length
            return True


    def print_targets(self):
        for ctg in self.targets:
            #print (ctg, self.cov_target)
            for target in self.targets[ctg]:
                localdepth = self.targets[ctg][target]['cov']/self.targets[ctg][target]['size']
                #print (target,self.targets[ctg][target]['cov'],self.targets[ctg][target]['size'],localdepth)
                if localdepth >= self.cov_target:
                    #print (self.conditions)
                    #print (target,ctg)
                    if ctg not in self.targets_to_skip.keys():
                        self.targets_to_skip[ctg]=set()
                    self.targets_to_skip[ctg].add(target)
                    print ("Try removing {}".format((self.targets[ctg][target]['cov'],self.targets[ctg][target]['cov']+self.targets[ctg][target]['size'])))
        print (self.targets_to_skip)

        target_list = list()
        for condition in self.conditions:
            for strand in condition.coords:
                for chr in condition.coords[strand]:
                    for coord in condition.coords[strand][chr]:
                        #print("{},{},{},{}".format(chr, coord[0], coord[1], strand))
                        if chr not in self.targets_to_skip.keys():
                            target_list.append("{},{},{},{}".format(chr, coord[0], coord[1], strand))
                        elif coord[0] not in self.targets_to_skip[chr]:
                            target_list.append("{},{},{},{}".format(chr, coord[0], coord[1], strand))

        self.write_new_toml(target_list, self.tomlfile)

    def refresh_index(self, reference):
        self.mappingobjects[reference]["last_used"] = time.time()
