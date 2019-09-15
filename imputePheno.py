import os
import sys
import glob
import pickle
import optparse
import pandas as pd
import numpy as np
import scipy.stats as st

class traits:
	traitFileName=""
	traitIndex = 0
	beta = 0
	sebeta = 0
	zScore = 0

	def __init__(self, name, index, beta, sebeta):
		self.traitFileName = name
		self.traitIndex = index
		self.beta = beta
		self.sebeta = sebeta
		self.zScore = float(beta)/float(sebeta)

def detectFileHander(input_folder, study_name, trait_names):
	file_handlers = []
	for trait in trait_names:
		file_name = input_folder + "/" + study_name + "_" + trait + ".txt"
		file_handler = open(file_name, "r")
		file_handlers.append(file_handler)
	return file_handlers

#compute the Z-score of phenotype to impute from all the other
#phenoypes and correlation. We assume the first phenotype (index 0)
#is the hard-to-collect phenotype that we want to impute.
#Assumption: First phenotype (index 0) to impute using ALL phenotype
def computePheno(R, Y):
	Sigma = R[1:,1:]
	Rl    = R[1:,0:1]
	Sigma_inv = np.linalg.inv(Sigma)
	return (np.matmul(np.matmul(Y, Sigma_inv),Rl))

# This function assume the first phenotypes is getting
# perdicted and all the other phenotypes are used to
# impute it.
def computeImpZAll (R, snp2z_map, output_file):
	output_file_hander = open(output_file, 'w')
	output_file_hander.write("SNP_ID\tZ_stat\tP_Val\n")
	for key,listTrait in snp2z_map.items():
		snpID = key
		ZList = []
		indexList = [0]
		for value in listTrait:
			ZList.append(value.zScore)
			indexList.append(value.traitIndex+1)
		Rtmp = R[indexList,:]
		Rtmp = Rtmp[:, indexList]
		if(len(ZList)==0):
			output_file_hander.write("%s\t%f\t%e\n" %(snpID, ZList[0] * Rtmp[0,1], p_value(ZList[0] * Rtmp[0,1])) )
		else:
			Z = computeImpZ(Rtmp, np.transpose(np.matrix(ZList))  )
			output_file_hander.write("%s\t%f\t%e\n"  %(snpID, Z[0,0], p_value(Z[0,0])) )
	output_file_hander.close()

def main (parser):
	(options, args) = parser.parse_args()
	cor_file = options.cor_file
	toimpute_name = options.toimpute_name
	input_file = options.input_file
	output_file = options.output_file
	trait_names = list(pd.read_csv(input_file, delim_whitespace=True).columns.values)[1:]
	print ("size %s" %(trait_names))	
	
	to_impute_index = -1
	all_phenotype_index = []
	phenotype_list = list(pd.read_csv(cor_file, delim_whitespace=True).columns.values)
	print (phenotype_list)
	R = np.loadtxt(open(cor_file), skiprows=1) # load the phenotype correlation between all phenotypes in the correlation file
	to_impute_index = phenotype_list.index(toimpute_name)
	all_phenotype_index.append(to_impute_index)
	for phen in trait_names:
		all_phenotype_index.append(phenotype_list.index(phen))	
	R = R[all_phenotype_index, :]
	R = R[:, all_phenotype_index]
	phen_data = np.loadtxt(open(input_file), skiprows=1)
	phen_mean = np.nanmean(phen_data[:,1:], axis=0)
	phen_var  = np.nanstd(phen_data[:,1:], axis=0)
	standadized_phen_data = (phen_data[:,1:]-phen_mean)/phen_var
	impY = computePheno (R, standadized_phen_data)
	outFileHandler = open(output_file, 'w')
	outFileHandler.write("SID %s %s\n" % (' '.join(trait_names[1:]), toimpute_name))
	for (ind_id, standadized_phen_val, imp_phen_val) in zip (phen_data[:,0], standadized_phen_data, impY):
		outFileHandler.write ("%d %s %f\n" %(ind_id, ' '.join(map(str, standadized_phen_val)), imp_phen_val))
	outFileHandler.close()

if __name__== "__main__":
	parser = optparse.OptionParser("usage: %prog [options] ")
	parser.add_option("-i", "--inputFile", dest="input_file",
    		default="", type="string",
			help="specify the File name that consist of all the phenotyes")
	parser.add_option("-f", "--phenotypeName", dest="toimpute_name",
            default="", type="string",
            help="specify the name of phenotype to impute")
	parser.add_option("-o", "--outputFile", dest="output_file",
			default="", type="string",
			help="specify the output File name")
	parser.add_option("-r", "--corFile", dest="cor_file",
			default="", type="string",
			help="specify the pair-wise correlation between different phenotype")
	main(parser)
