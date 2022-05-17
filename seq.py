import os
from Bio import SeqIO
from kmer import kmer_featurization
import numpy as np
import sys


if len(sys.argv) == 3:
    if sys.argv[1] == "--k":
        if int(sys.argv[2]) in [1,2,3,4,5,6,7,8,9]:
            k = int(sys.argv[2])
            fastaData = []
            for filename in os.listdir(os.getcwd() + "/fastas/"):
                intermediate = SeqIO.parse(os.getcwd() + "/fastas/" + filename,'fasta')
                for item in intermediate:
                    fastaData.append(item)
            obj = kmer_featurization(k)
            kmer_features = obj.obtain_kmer_feature_for_a_list_of_sequences(fastaData,write_number_of_occurrences=False)        
            print(kmer_features)
else:
    print("bitte per --k anzahl k-mer angeben")