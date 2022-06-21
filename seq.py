import os
from Bio import SeqIO
from Bio.Cluster import pca
from kmer import kmer_featurization
import numpy as np
import sys
import matplotlib.pyplot as plt

if len(sys.argv) == 3:
	if sys.argv[1] == "-k":
		if int(sys.argv[2]) in [1,2,3,4,5,6,7,8,9]:
			k = int(sys.argv[2])
			fastaData = []
			matrixpca = [[]]
			i = 0
			for filename in os.listdir(os.getcwd() + "/fastas/"):
				intermediate = SeqIO.parse(os.getcwd() + "/fastas/" + filename,'fasta')
				j = 0
				for item in intermediate:
					fastaData.append(item)
			obj = kmer_featurization(k)
			kmer_features = obj.obtain_kmer_feature_for_a_list_of_sequences(fastaData,write_number_of_occurrences=False)      
			matrixpca = kmer_features
			columnmean, coordinates, components, eigenvalues = pca(matrixpca)
			xco = []
			yco = []
			i = 0
			for filename in os.listdir(os.getcwd() + "/fastas/"):
				xco.append(coordinates[i][0])
				yco.append(coordinates[i][1])
				i = i + 1
			plt.scatter(xco, yco)
			#erster eintrage von jedem array x und zweiter einrag von jedem array  y nicht 
			plt.show()
 
else:
    print("bitte per -k anzahl k-mer angeben")

