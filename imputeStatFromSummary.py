import os
import sys
import glob
import pickle
import optparse
import pandas as pd
import numpy as np
import scipy.stats as st

class traits:
	traitFileName="";
	traitIndex = 0;
	beta = 0;
	sebeta = 0;
	zScore = 0;

	def __init__(self, name, index, beta, sebeta):
		self.traitFileName = name;
		self.traitIndex = index;
		self.beta = beta;
		self.sebeta = sebeta;
		self.zScore = float(beta)/float(sebeta);

## This three functions are used to save or load
## python object as pickle files
### BEGIN Picle utilities ####
def initalize_pickle_folder(name):
	if not os.path.exists(name):
		path = os.getcwd()
		try:
			os.mkdir(path + "/obj")
		except OSError:
			print ("Creation of the directory %s failed" % path)		

def save_obj(obj, name):
	with open('obj/'+ name + '.pkl', 'wb') as f:
		pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
	print ("Save data as Pickle file");

def load_obj(name):
	print ("Load data from Pickle file");	
	with open('obj/' + name + '.pkl', 'rb') as f:
		return pickle.load(f)
### FINISH Picle utilities ####

def detectFileHander(input_folder, study_name, trait_names):
	file_handlers = []
	for trait in trait_names:
		file_name = input_folder + "/" + study_name + "_" + trait + ".txt"
		file_handler = open(file_name, "r")
		file_handlers.append(file_handler)
	return file_handlers

def getSNP2ZDictionary(file_handlers, trait_names):
	snp2z_map = {};
	index = 0;
	for handler in file_handlers:
		for line in handler:
			if "SNP" in line:
				continue;
			data = line.split("\t");
			snpID  = data[0];
			beta   = data[7];
			sebeta = data[8];
			if snpID in snp2z_map:
				snp2z_map[snpID].append(traits(trait_names[index], index, beta, sebeta));
			else:
				tmpList = [];
				tmpList.append(traits(trait_names[index], index, beta, sebeta));
				snp2z_map[snpID] = tmpList;
		index = index + 1;
	return snp2z_map; 

def p_value(stat) :
	return(2 * st.norm.pdf(-abs(stat), loc=0, scale=1));

#compute the Z-score of phenotype to impute from all the other
#phenoypes and correlation. We assume the first phenotype (index 0)
#is the hard-to-collect phenotype that we want to impute.
#Assumption: First phenotype (index 0) to impute using ALL phenotype
def computeImpZ(R, Z):
	Sigma = R[1:,1:];
	Rl    = R[1:,0:1]
	return (np.dot(np.transpose(np.dot(np.linalg.inv(Sigma), Rl)), Z));

# This function assume the first phenotypes is getting
# perdicted and all the other phenotypes are used to
# impute it.
def computeImpZAll (R, snp2z_map, output_file):
	output_file_hander = open(output_file, 'w')
	output_file_hander.write("SNP_ID\tZ_stat\tP_Val\n")
	for key,listTrait in snp2z_map.items():
		snpID = key;
		ZList = [];
		indexList = [0];
		for value in listTrait:
			ZList.append(value.zScore);
			indexList.append(value.traitIndex+1);
		Rtmp = R[indexList,:]
		Rtmp = Rtmp[:, indexList];
		if(len(ZList)==0):
			output_file_hander.write("%s\t%f\t%e\n" %(snpID, ZList[0] * Rtmp[0,1], p_value(ZList[0] * Rtmp[0,1])) )
		else:
			Z = computeImpZ(Rtmp, np.transpose(np.matrix(ZList))  );
			output_file_hander.write("%s\t%f\t%e\n"  %(snpID, Z[0,0], p_value(Z[0,0])) )
	output_file_hander.close()

def check_files(phenotype_list, trait_names, input_folder):
	for name in trait_names:
		if not name in phenotype_list:
			print (f"We expect a file for {name} in your {input_folder}")
			raise Exception("At least one phenotype file is missing")

def main (parser):
	(options, args) = parser.parse_args()
	cor_file      =  options.cor_file
	study_name    =  options.study_name
	toimpute_name =  options.toimpute_name
	input_folder  = options.input_folder
	output_file   = options.output_file

	trait_names = [w.replace(input_folder+ "/" + study_name + "_", '').split(".")[0] for w in glob.glob(input_folder + "/" + study_name + "*")]

	if (len(trait_names) <= 1):
		print (f"The {input_folder} does not have enough phenotypes or you given the wrong path")
		sys.exit(0)
	
	snp2z_map   = {}
	to_impute_index = -1
	all_phenotype_index = []

	file_handlers = detectFileHander(input_folder, study_name, trait_names) 
	initalize_pickle_folder(study_name)
	if not (os.path.exists('obj/'+study_name + ".pkl")):
		snp2z_map = getSNP2ZDictionary(file_handlers, trait_names);
		save_obj(snp2z_map, study_name);
	else:
		snp2z_map = load_obj(study_name);
	print ("size %i" %(len(snp2z_map)))	
	
	phenotype_list = list(pd.read_csv(cor_file, delim_whitespace=True).columns.values)
	try:
		check_files(phenotype_list, trait_names, input_folder)
	except Exception as error:
		print (err)	
		sys.exit(0)		
	
	R = np.loadtxt(open(cor_file), skiprows=1) # load the phenotype correlation between all phenotypes in the correlation file
	
	to_impute_index = phenotype_list.index(toimpute_name)
	all_phenotype_index.append(to_impute_index);
	for phen in trait_names:
		all_phenotype_index.append(phenotype_list.index(phen))	

	R = R[all_phenotype_index, :];
	R = R[:, all_phenotype_index];

	computeImpZAll(R, snp2z_map, output_file);

if __name__== "__main__":
	parser = optparse.OptionParser("usage: %prog [options] ")
	parser.add_option("-i", "--inputFolder", dest="input_folder",
    		default="", type="string",
			help="specify the Folder name that consist of all the phenotyes")
	parser.add_option("-n", "--studyName", dest="study_name",
            default="", type="string",
            help="specify the name of GWAS study")
	parser.add_option("-f", "--phenotypeName", dest="toimpute_name",
            default="", type="string",
            help="specify the name of phenotype to impute")
	parser.add_option("-o", "--outputFile", dest="output_file",
			default="", type="string",
			help="specify the output File name");
	parser.add_option("-r", "--corFile", dest="cor_file",
			default="", type="string",
			help="specify the pair-wise correlation between different phenotype")
	main(parser)
