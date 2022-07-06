import os
from Bio import SeqIO
from Bio.Cluster import pca
from kmer import kmer_featurization
import numpy as np
import sys
import matplotlib.pyplot as plt
import argparse
xco = []
yco = []
fastaData = []
matrixpca = [[]]
fasta_names = []
i = 0
j = 0
k = 0
parser = argparse.ArgumentParser("This programm calculates the k-mers of choosen Fasta Data and plots a 2d pca. FASTA Data has to be saved to fasta directory and number of k-mers has to be given as command argument -k.")
parser.add_argument("-k","--k_mers", help="number of k-mers", type=int, choices=[1,2,3,4,5,6,7,8,9])
args = parser.parse_args()
    
def k_mer_calc():
	obj = kmer_featurization(args.k_mers)
	return obj.obtain_kmer_feature_for_a_list_of_sequences(fastaData,write_number_of_occurrences=False)
    	
def plot_pca(matrixpca):
	columnmean, coordinates, components, eigenvalues = pca(matrixpca)
	i = 0
	for filename in os.listdir(os.getcwd() + "/fastas/"):
		xco.append(coordinates[i][0])
		yco.append(coordinates[i][1])
		i = i + 1
	plt.scatter(xco, yco)
	i = 0
	for i, txt in enumerate(fasta_names):
		plt.annotate(txt, (xco[i],yco[i]))
	plt.show()    	
    
def main():
	k = int(args.k_mers)
	for filename in os.listdir(os.getcwd() + "/fastas/"):
		fasta_names.append(filename)
		intermediate = SeqIO.parse(os.getcwd() + "/fastas/" + filename,'fasta')
		for item in intermediate:
			fastaData.append(item)				
	plot_pca(k_mer_calc())
	
if __name__ == "__main__":
    main()	
 


