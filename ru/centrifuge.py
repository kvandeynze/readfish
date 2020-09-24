"""
A helper function to provide a class that can run centrifuge and return results as required by iteralign
"""
import subprocess
import csv


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

    def add_taxon(self,taxID,name,genomeSize,numUniqueReads):
        if taxID not in self.tax_data.keys():
            self.tax_data[taxID]=dict()
            self.tax_data[taxID]["name"]=name
            self.tax_data[taxID]["genomeSize"]=genomeSize
            self.tax_data[taxID]["uniqueReads"]=int(numUniqueReads)
        else:
            self.tax_data[taxID]["uniqueReads"]+=int(numUniqueReads)

    def _calculate_targets(self):
        for taxID in self.tax_data:
            if self.tax_data[taxID]["uniqueReads"]>=self.threshold:
                self.new_targets.append(taxID)


    def classify(self,fastqfileList):
        # convert the 'fastqfileList' into a string valid for the list of fastq files to be read by centrifuge
        fastq_str = ",".join(fastqfileList)

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

        #First grab and store information on the genomes seen.
        with open(self.args.creport, newline='') as f:
            reader = csv.reader(f)
            next(reader)
            for row in reader:
                #print(row)
                name, taxID, taxRank, genomeSize, numReads, numUniqueReads, abundance = row[0].split("\t")
                self.add_taxon(taxID,name,genomeSize,numUniqueReads)

        """
        #This bit of code is only needed if we need access to the readIDs - which we do not!
        for line in out.splitlines()[1:-1]:
            #print (line)
            readID, seqID, taxID, score , secondBestScore, hitLength,queryLength, numMatches= line.split("\t")
            if numMatches == 1: #this filters out reads that map to one or more genomes
                print (readID,numMatches,taxID)
        """
        print (self.tax_data)
        self._calculate_targets()
        ### So if we have new targets we need to download the references and get them for building an index.



        #name,taxID,taxRank,genomeSize,numReads,numUniqueReads,abundance =