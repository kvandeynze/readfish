"""
A helper function to provide a class that can run centrifuge and return results as required by iteralign
"""
import subprocess
import csv
import urllib.request as request
import urllib.error as url_error
from io import StringIO, BytesIO
import gzip
from Bio import SeqIO
import logging
import logging.handlers



class CentrifugeServer():
    """
    Long term goal is that this will allow streaming of individual reads or chunks of reads to a continuously open centrifuge classifier.
    """

    def __init__(self,args):
        print ("Configuring centrifuge")
        self.args = args
        self.tax_data = dict()
        self.threshold = self.args.threshold
        self.new_targets = list()
        self.all_targets = set()
        self.ref_lookup = dict()
        self._store_urls()
        self.plasmidsd = dict()
        if self.args.plasmids:
            self.plasmidsd = self._store_plasmids()
        # ToDO Do we need this? Perhaps - it accepts a list of sequences that you might want to reject up front.
        if self.args.references:
            for taxID in self.args.references:
                self.new_targets.append(str(taxID))


    def _store_plasmids(self):
        r = ("name", "path")
        with open(self.args.csummary) as f:
            d = {int(x[0]): dict(zip(r, x[1:])) for i, l in enumerate(f) for x in (l.strip().split("\t"),)
                 if i > 0}

        return d


    def _store_urls(self):
        with open(self.args.csummary, newline='') as f:
            reader = csv.reader(f)
            next(reader)
            for row in reader:
                taxID, name,url = row[0].split("\t")
                #print (taxID,name,url)
                self.ref_lookup[taxID]=url


    def _download_references(self, taxID):
        lengths = {}
        link = self.ref_lookup[taxID]
        print("Attempting to download: {}".format(link), end="\n")
        try:
            response = request.urlopen(link)
        except url_error.URLError as e:
            print(e)
            print("Closing script")

        compressed_file = BytesIO(response.read())
        with gzip.open(compressed_file, "rt") as fh, gzip.open(self.args.path + self.args.prefix + self.args.gfasta,
                                                               "at") as fasta_seq:
           for seq_record in SeqIO.parse(fh, "fasta"):
                if len(seq_record) > self.args.seqlength:
                    lengths[seq_record.id] = len(seq_record)
                    SeqIO.write(seq_record, fasta_seq, "fasta")

        if self.args.plasmids:
            logging.info("Obtaining the plasmids for the following taxids: {}".format(taxID))

            with gzip.open(self.args.plasmids, "rt") as fh, gzip.open(self.args.path + self.args.prefix + self.args.gfasta,
                                                                         "at") as fasta_seq:
                for seq_record in SeqIO.parse(fh, "fasta"):
                    if taxID in self.plasmidsd.keys():
                        if self.plasmidsd[taxID]["name"] in seq_record.description and len(seq_record) > self.args.seqlength:
                            lengths[seq_record.id] = len(seq_record)
                            SeqIO.write(seq_record, fasta_seq, "fasta")

        logging.info("Genome file generated in {}".format(fasta_seq))

        return lengths

    def _add_taxon(self,taxID,name,genomeSize,numUniqueReads):
        if taxID not in self.tax_data.keys():
            self.tax_data[taxID]=dict()
            self.tax_data[taxID]["name"]=name
            self.tax_data[taxID]["genomeSize"]=genomeSize
            self.tax_data[taxID]["uniqueReads"]=int(numUniqueReads)
        else:
            self.tax_data[taxID]["uniqueReads"]+=int(numUniqueReads)

    def _calculate_targets(self):
        for taxID in self.tax_data:
            if taxID not in self.all_targets and self.tax_data[taxID]["uniqueReads"]>=self.threshold:
                self.new_targets.append(taxID)
                self.all_targets.add(taxID)

    def classify(self,fastqfileList,mapper):
        # convert the 'fastqfileList' into a string valid for the list of fastq files to be read by centrifuge
        fastq_str = ",".join(fastqfileList)

        print ("Reference is {}.".format(self.args.toml['conditions']['reference']))

        # centrifuge command to classify reads in the fastq files found by watchdog
        centrifuge_cmd = "centrifuge -p {} -x {} -q {}".format(self.args.threads,
                                                               self.args.cindex,
                                                               fastq_str
                                                               )


        # subprocess for 'centrifuge_cmd'
        proc = subprocess.Popen(
            centrifuge_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            stdin=subprocess.PIPE,
            shell=True,
            # Aliased by `text=True` in 3.7
            universal_newlines=True,
        )
        out, err = proc.communicate()
        proc.stdout.close()

        print ("finished that subprocess")

        #First grab and store information on the genomes seen.
        with open(self.args.creport, newline='') as f:
            reader = csv.reader(f)
            next(reader)
            for row in reader:
                #print(row)
                name, taxID, taxRank, genomeSize, numReads, numUniqueReads, abundance = row[0].split("\t")
                self._add_taxon(taxID,name,genomeSize,numUniqueReads)

        """
        #This bit of code is only needed if we need access to the readIDs - which we do not!
        for line in out.splitlines()[1:-1]:
            #print (line)
            readID, seqID, taxID, score , secondBestScore, hitLength,queryLength, numMatches= line.split("\t")
            if numMatches == 1: #this filters out reads that map to one or more genomes
                print (readID,numMatches,taxID)
        """
        #print (self.tax_data)
        self._calculate_targets()

        ### So if we have new targets we need to download the references and get them for building an index.
        if self.new_targets:
            while len(self.new_targets) > 0:
                target = self.new_targets.pop()
                print(self._download_references(target))

            #Make a new reference

            mapper.load_index("test",self.args.path + self.args.prefix + self.args.gfasta)
            self.args.toml['conditions']['reference']=self.args.path + self.args.prefix + self.args.gfasta
            print ("Updated reference is {}".format(self.args.toml['conditions']['reference']))
        #return file handle and instruct to make a new reference IF we need to do that.




        #name,taxID,taxRank,genomeSize,numReads,numUniqueReads,abundance =