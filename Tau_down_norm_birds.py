#Filter rMATS output JC files

#!/Users/thearogers/anaconda2/bin/python
# -*- coding: utf-8 -*- 
#==============================================================================
import argparse
import sys
from collections import defaultdict
import os, time
import numpy as np 
import math
import re
import random
from scipy.stats import wilcoxon
#==============================================================================
#Command line options==========================================================
#==============================================================================

#============================
parser = argparse.ArgumentParser()

parser.add_argument("SE_files", type=str,nargs=1, help="SE.MATS.JC.txt file") 

parser.add_argument("MXE_files", type=str,nargs=1, help="MXE.MATS.JC.txt file") 

parser.add_argument("readcoverage", type=str, help="Number of reads required to support splice junction")

parser.add_argument("coordinates", type=str,
					help="Tab-separated file containing ortholog sex-bias info; for example \'ENSAPLG00000021936	MB\'")

if len(sys.argv)==1:
	parser.print_help()
	sys.exit(1)
args = parser.parse_args()

#==============================================================================
#Functions
#==============================================================================

def filter_SE(name, readcoverage):
	with open(name, "r") as rmats_files:
		gene_list =[]
		pass_filt = []
		a = []
		b = []
		c = []
		d = []
		for line in rmats_files:
			if not line.startswith("ID"):
				gene_id = line.split("\t")[1].split('"')[1]
				gene_list.append(gene_id)

				ijc_sample1 = line.split("\t")[12].split(",")
				ijc_sample1 = [float(i) for i in ijc_sample1]
				ijc_sample1_threshold = len([i for i in ijc_sample1 if float(i)>=(float(args.readcoverage)/2)]) 

				sjc_sample1 = line.split("\t")[13].split(",")
				sjc_sample1 = [float(i) for i in sjc_sample1]
				sjc_sample1_threshold = len([i for i in sjc_sample1 if float(i)>=(float(args.readcoverage)/2)])

				sample1_rc_sum = (np.add(ijc_sample1,sjc_sample1).tolist())
				sample1_threshold = len([i for i in sample1_rc_sum if i>=(float(args.readcoverage))]) 

				ijc_sample2 = line.split("\t")[14].split(",")
				ijc_sample2 = [float(i) for i in ijc_sample2]
				ijc_sample2_threshold = len([i for i in ijc_sample2 if float(i)>=(float(args.readcoverage)/2)]) 

				sjc_sample2 = line.split("\t")[15].split(",")
				sjc_sample2 = [float(i) for i in sjc_sample2]
				sjc_sample2_threshold = len([i for i in sjc_sample2 if float(i)>=(float(args.readcoverage)/2)]) 

				sample2_rc_sum = (np.add(ijc_sample2,sjc_sample2).tolist())
				sample2_threshold = len([i for i in sample2_rc_sum if i>=(float(args.readcoverage))]) 
				

				if (sample1_threshold > 2) and (sample2_threshold > 2):
					if ((ijc_sample1_threshold > len(ijc_sample1)/2) and (sjc_sample1_threshold > len(sjc_sample1)/2)) or ((ijc_sample1_threshold > len(ijc_sample1)/2) and (sjc_sample2_threshold > len(sjc_sample2)/2)) or ((sjc_sample1_threshold > len(sjc_sample1)/2) and (ijc_sample2_threshold > len(ijc_sample2)/2)) or (ijc_sample2_threshold > len(ijc_sample2)/2) and (sjc_sample2_threshold > len(sjc_sample2)/2):

						pass_filt.append(gene_id)
						a.append(ijc_sample1)
						b.append(sjc_sample1)
						c.append(ijc_sample2)
						d.append(sjc_sample2)



		print "Number of SE genes = ", len(set(gene_list))
		print "Number of SE junctions = ", len(gene_list)

		return pass_filt, a, b, c, d 
					
			
def filter_MXE(name,readcoverage):
	with open(name, "r") as rmats_files:
		gene_list =[]
		pass_filt = []
		a = []
		b = []
		c = []
		d = []
		i=0

		for line in rmats_files:
			if not line.startswith("ID"): 
				gene_id = line.split("\t")[1].split('"')[1]
				gene_list.append(gene_id)

				ijc_sample1 = line.split("\t")[14].split(",")
				ijc_sample1 = [float(i) for i in ijc_sample1]
				ijc_sample1_threshold = len([i for i in ijc_sample1 if float(i)>=(float(args.readcoverage)/2)]) 

				sjc_sample1 = line.split("\t")[15].split(",")
				sjc_sample1 = [float(i) for i in sjc_sample1]
				sjc_sample1_threshold = len([i for i in sjc_sample1 if float(i)>=(float(args.readcoverage)/2)]) 

				sample1_rc_sum = (np.add(ijc_sample1,sjc_sample1).tolist())
				sample1_threshold = len([i for i in sample1_rc_sum if i>=(float(args.readcoverage))]) 

				ijc_sample2 = line.split("\t")[16].split(",")
				ijc_sample2 = [float(i) for i in ijc_sample2]
				ijc_sample2_threshold = len([i for i in ijc_sample2 if float(i)>=(float(args.readcoverage)/2)]) 

				sjc_sample2 = line.split("\t")[17].split(",")
				sjc_sample2 = [float(i) for i in sjc_sample2]
				sjc_sample2_threshold = len([i for i in sjc_sample2 if float(i)>=(float(args.readcoverage)/2)]) 

				sample2_rc_sum = (np.add(ijc_sample2,sjc_sample2).tolist())
				sample2_threshold = len([i for i in sample2_rc_sum if i>=(float(args.readcoverage))]) 

				if (sample1_threshold > 2) and (sample2_threshold > 2):	
					if ((ijc_sample1_threshold > len(ijc_sample1)/2) and (sjc_sample1_threshold > len(sjc_sample1)/2)) or ((ijc_sample1_threshold > len(ijc_sample1)/2) and (sjc_sample2_threshold > len(sjc_sample2)/2)) or ((sjc_sample1_threshold > len(sjc_sample1)/2) and (ijc_sample2_threshold > len(ijc_sample2)/2)) or (ijc_sample2_threshold > len(ijc_sample2)/2) and (sjc_sample2_threshold > len(sjc_sample2)/2):
						pass_filt.append(gene_id)
						a.append(ijc_sample1)
						b.append(sjc_sample1)
						c.append(ijc_sample2)
						d.append(sjc_sample2)
			
		print "Number of MXE genes = ", len(set(gene_list))
		print "Number of MXE junctions = ", len(gene_list)

		return pass_filt, a, b, c, d 



def get_lib_size(passfilt, a, b, c, d):

	s1_1_ijc = []
	s1_1_sjc = []
	s1_2_ijc = []
	s1_2_sjc = []
	s1_3_ijc = []
	s1_3_sjc = []
	s1_4_ijc = []
	s1_4_sjc = []

	s2_1_ijc = []
	s2_1_sjc = []
	s2_2_ijc = []
	s2_2_sjc = []
	s2_3_ijc = []
	s2_3_sjc = []
	s2_4_ijc = []
	s2_4_sjc = []

	#i = 0

	for i in range(0,len(passfilt)):
		#print(a[i])

		s1_1_ijc.append(int(float(a[i][0])))
		s1_2_ijc.append(int(float(a[i][1])))
		s1_3_ijc.append(int(float(a[i][2])))
		s1_4_ijc.append(int(float(a[i][3])))

		s1_1_sjc.append(int(float(b[i][0])))
		s1_2_sjc.append(int(float(b[i][1])))
		s1_3_sjc.append(int(float(b[i][2])))
		s1_4_sjc.append(int(float(b[i][3])))

		s2_1_ijc.append(int(float(c[i][0])))
		s2_2_ijc.append(int(float(c[i][1])))
		s2_3_ijc.append(int(float(c[i][2])))
		s2_4_ijc.append(int(float(c[i][3])))

		s2_1_sjc.append(int(float(d[i][0])))
		s2_2_sjc.append(int(float(d[i][1])))
		s2_3_sjc.append(int(float(d[i][2])))
		s2_4_sjc.append(int(float(d[i][3])))

		s1_1_sum = sum(s1_1_ijc)+sum(s1_1_sjc)
		s1_2_sum = sum(s1_2_ijc)+sum(s1_2_sjc)
		s1_3_sum = sum(s1_3_ijc)+sum(s1_3_sjc)
		s1_4_sum = sum(s1_4_ijc)+sum(s1_4_sjc)

		s2_1_sum = sum(s2_1_ijc)+sum(s2_1_sjc)
		s2_2_sum = sum(s2_2_ijc)+sum(s2_2_sjc)
		s2_3_sum = sum(s2_3_ijc)+sum(s2_3_sjc)
		s2_4_sum = sum(s2_4_ijc)+sum(s2_4_sjc)

		#print s1_1_sum, s1_2_sum, s1_3_sum, s1_4_sum, s2_1_sum, s2_2_sum, s2_3_sum, s2_4_sum



	return s1_1_sum, s1_2_sum, s1_3_sum, s1_4_sum, s2_1_sum, s2_2_sum, s2_3_sum, s2_4_sum

def get_gene_max(genes, a, b, c, d):	

	gene_max_s1 = {}
	gene_max_s2 = {}

	for i in range(0,len(genes)):
		gene = genes[i]
		max_count= max(sum(a[i]), sum(b[i]))
		if gene in gene_max_s1:
			gene_max_s1[gene] = max(max_count,gene_max_s1[gene])
		else:
			gene_max_s1[gene] = max_count
	for i in range(0,len(genes)):
		gene = genes[i]
		max_count= max(sum(c[i]), sum(d[i]))
		if gene in gene_max_s2:
			gene_max_s2[gene] = max(max_count,gene_max_s2[gene])
		else:
			gene_max_s2[gene] = max_count
			
	return gene_max_s1, gene_max_s2

def get_tau(genes, a, b, c, d, gm_s1, gm_s2, gc_s1, gc_s2):
	tau_dict_s1 = {}
	tau_dict_s2 = {}
	for i in range(0,len(genes)):
		gene = genes[i]
		ijc_s1 = sum(a[i])
		sjc_s1 = sum(b[i])
		ijc_s2 = sum(c[i])
		sjc_s2 = sum(d[i])
		#print "ijcs1", ijc_s1

		#print gene, a[gene], b[gene], gm_s1[gene], c[gene], d[gene], gm_s2[gene]
		x_i_s1 = (1-(math.log(float(ijc_s1+1.00))/math.log(float(gm_s1[gene]+1.00)))) + (1-(math.log(float(sjc_s1+1.00))/math.log(float(gm_s1[gene]+1.00))))
		x_i_s2 = (1-(math.log(float(ijc_s2+1.00))/math.log(float(gm_s2[gene]+1.00)))) + (1-(math.log(float(sjc_s2+1.00))/math.log(float(gm_s2[gene]+1.00))))
		#if gene == "FBgn0031187":
			#print gene, ijc_s1, sjc_s1, gm_s1[gene], ijc_s2, sjc_s2, gm_s2[gene]	
		if gene in tau_dict_s1:
			tau_dict_s1[gene] += x_i_s1
		else:
			tau_dict_s1[gene] = x_i_s1

		if gene in tau_dict_s2:
			tau_dict_s2[gene] += x_i_s2
		else:
			tau_dict_s2[gene] = x_i_s2


	tau_s1 = {}
	tau_s2 = {}
	
	for gene in tau_dict_s1:
	 	raw_tau_s1 = tau_dict_s1[gene]
	 	isoform_no_s1 = len(gc_s1[gene]) - 1
	 	tau_s1[gene] = (raw_tau_s1/isoform_no_s1)

	 	raw_tau_s2 = tau_dict_s2[gene]
	 	isoform_no_s2 = len(gc_s2[gene]) - 1
	 	tau_s2[gene] = (raw_tau_s2/isoform_no_s2)
	 	
	return tau_s1, tau_s2


def downsample(passfilt, a, b, c, d):
	a_down = []
	b_down = []
	c_down = []
	d_down = []
	

	for i in range(0,len(passfilt)):
		gene_id = passfilt[i]
		
		if (sum(a[i])+sum(b[i])) > (sum(c[i])+sum(d[i])):
			m_f = (sum(a[i])+sum(b[i]))/(sum(c[i])+sum(d[i]))
			l_a_down = [(x/m_f) for x in a[i]]
			l_b_down = [(x/m_f) for x in b[i]]
			a_down.append(l_a_down)
			b_down.append(l_b_down)
			c_down.append(c[i])
			d_down.append(d[i])
			
		elif (sum(c[i])+sum(d[i])) > (sum(a[i])+sum(b[i])):
			f_m = (sum(c[i])+sum(d[i]))/(sum(a[i])+sum(b[i]))
			l_c_down = [(x/f_m) for x in c[i]]
			l_d_down = [(x/f_m) for x in d[i]]
			c_down.append(l_c_down)
			d_down.append(l_d_down)
			a_down.append(a[i])
			b_down.append(b[i])

		else:
			a_down.append(a[i])
			b_down.append(b[i])
			c_down.append(c[i])
			d_down.append(d[i])

	return a_down, b_down, c_down, d_down

def normalize(passfilt, a, b, c, d):
	a_norm = []
	b_norm = []
	c_norm = []
	d_norm = []
	a_norm_sum = []
	b_norm_sum = []
	c_norm_sum = []
	d_norm_sum = []

	s1_1_lib, s1_2_lib, s1_3_lib, s1_4_lib, s2_1_lib, s2_2_lib, s2_3_lib, s2_4_lib = get_lib_size(passfilt, a, b, c, d)

	#print s1_1_lib, s1_2_lib, s1_3_lib, s1_4_lib, s2_1_lib, s2_2_lib, s2_3_lib, s2_4_lib

	s1_1_lib_norm = 0
	s1_2_lib_norm = 0

	for i in range(0, len(passfilt)):
		l_a = []
		l_b = []
		l_c = []
		l_d = []
		
		l_a = ((float(a[i][0])/float(s1_1_lib)), (float(a[i][1])/float(s1_2_lib)), (float(a[i][2])/float(s1_3_lib)), (float(a[i][3])/float(s1_4_lib)))
		l_b = ((float(b[i][0])/float(s1_1_lib)), (float(b[i][1])/float(s1_2_lib)), (float(b[i][2])/float(s1_3_lib)), (float(b[i][3])/float(s1_4_lib)))
		l_c = ((float(c[i][0])/float(s2_1_lib)), (float(c[i][1])/float(s2_2_lib)), (float(c[i][2])/float(s2_3_lib)), (float(c[i][3])/float(s2_4_lib)))
		l_d = ((float(d[i][0])/float(s2_1_lib)), (float(d[i][1])/float(s2_2_lib)), (float(d[i][2])/float(s2_3_lib)), (float(d[i][3])/float(s2_4_lib)))
		
		a_norm.append(l_a)
		a_norm_sum.append(sum([x for x in l_a]))

		b_norm.append(l_b)
		b_norm_sum.append(sum([x for x in l_b]))

		c_norm.append(l_c)
		c_norm_sum.append(sum([x for x in l_c]))

		d_norm.append(l_d)
		d_norm_sum.append(sum([x for x in l_d]))

		s1_1_lib_norm += (l_a[0] + l_b[0])
		s1_2_lib_norm += (l_a[1] + l_b[1])


	ijc_s1_norm = {}
	sjc_s1_norm = {}
	ijc_s2_norm = {}
	sjc_s2_norm = {}

	for i in range(0,len(passfilt)):
		gene_id = passfilt[i]
		if gene_id in ijc_s1_norm:
			ijc_s1_norm[i] = np.add(ijc_s1_norm[i],a_norm[i])
		else:
			ijc_s1_norm[i] = a_norm[i]

		if gene_id in sjc_s1_norm:
			sjc_s1_norm[i] = np.add(sjc_s1_norm[i],b_norm[i])
		else:
			sjc_s1_norm[i] = b_norm[i]

		if gene_id in ijc_s2_norm:
			ijc_s2_norm[i] = np.add(ijc_s2_norm[i],c_norm[i])
		else:
			ijc_s2_norm[i] = c_norm[i]

		if gene_id in sjc_s2_norm:
			sjc_s2_norm[i] = np.add(sjc_s2_norm[i],d_norm[i])
		else:
			sjc_s2_norm[i] = d_norm[i]

	return ijc_s1_norm, sjc_s1_norm, ijc_s2_norm, sjc_s2_norm

def make_dictionary(passfilt, a, b, c, d):

	a_dict = {}
	b_dict = {}
	c_dict = {}
	d_dict = {}

	for i in range(0,len(passfilt)):
		gene_id = passfilt[i]
		if gene_id in a_dict:
			a_dict[gene_id] = a_dict[gene_id] + sum(a[i])
		else:
			a_dict[gene_id] = sum(a[i])

		if gene_id in b_dict:
			b_dict[gene_id] = b_dict[gene_id] + sum(b[i])
		else:
			b_dict[gene_id] = sum(b[i])

		if gene_id in c_dict:
			c_dict[gene_id] = c_dict[gene_id] + sum(c[i])
		else:
			c_dict[gene_id] = sum(c[i])

		if gene_id in d_dict:
			d_dict[gene_id] = d_dict[gene_id] + sum(d[i])
		else:
			d_dict[gene_id] = sum(d[i])

	return a_dict, b_dict, c_dict, d_dict

def bootstrap(list):
	length = 50
	tau = []
	bootstrap = 0
	while bootstrap < 1000:
		bootstrap += 1
		bootstrap_list = []
		for i in range(length):
			bootstrap_list.append(random.choice(list))
		if len(bootstrap_list) != length:
			print "ERROR - incorrect number sampled"
		boot_tau = np.mean(bootstrap_list)
		tau.append(float(boot_tau))
	

	lower_tau, upper_tau = get_confidence_intervals(tau)
	return lower_tau, upper_tau

def get_confidence_intervals(list):
    upper = round((np.percentile(list,97.5)),3)
    lower = round((np.percentile(list,2.5)),3)
    return lower, upper	

def calculate_stats(list, name):
	print "Calculating stats for ", name
	tau = np.mean(list)
	lower_tau, upper_tau = bootstrap(list)
	lower_upper_tau = str(round((lower_tau),3))+"-"+str(round((upper_tau),3))
	print name, len(list), round((tau),3), lower_upper_tau
	if name == "Male-biased expression s1":
		MB_s1_string = "ggplot()+geom_bar(aes(x=1,y="+str(round((tau),3))+",alpha=0.5),stat=\"identity\",fill=\"dodgerblue4\")+geom_errorbar(aes(x=1,ymin="+str(round((lower_tau),3))+",ymax="+str(round((upper_tau),3))+",width=.2),color=\"dodgerblue4\")+"
		return MB_s1_string
	elif name == "Male-biased expression s2":
		MB_s2_string = "geom_bar(aes(x=2,y="+str(round((tau),3))+",alpha=0.5),stat=\"identity\",fill=\"firebrick3\")+geom_errorbar(aes(x=2,ymin="+str(round((lower_tau),3))+",ymax="+str(round((upper_tau),3))+",width=.2),color=\"firebrick3\")+"
		return MB_s2_string
	elif name == "Female-biased expression s1":
		FB_s1_string = "geom_bar(aes(x=3,y="+str(round((tau),3))+",alpha=0.5),stat=\"identity\",fill=\"dodgerblue4\")+geom_errorbar(aes(x=3,ymin="+str(round((lower_tau),3))+",ymax="+str(round((upper_tau),3))+",width=.2),color=\"dodgerblue4\")+"
		return FB_s1_string
	elif name == "Female-biased expression s2":
		FB_s2_string = "geom_bar(aes(x=4,y="+str(round((tau),3))+",alpha=0.5),stat=\"identity\",fill=\"firebrick3\")+geom_errorbar(aes(x=4,ymin="+str(round((lower_tau),3))+",ymax="+str(round((upper_tau),3))+",width=.2),color=\"firebrick3\")+"
		return FB_s2_string
	elif name == "Unbiased expression s1":
		UB_s1_string = "geom_bar(aes(x=5,y="+str(round((tau),3))+",alpha=0.5),stat=\"identity\",fill=\"grey\")+geom_errorbar(aes(x=5,ymin="+str(round((lower_tau),3))+",ymax="+str(round((upper_tau),3))+",width=.2),color=\"grey\")+"
		return UB_s1_string
	elif name == "Unbiased expression s2":
		UB_s2_string = "geom_bar(aes(x=6,y="+str(round((tau),3))+",alpha=0.5),stat=\"identity\",fill=\"grey\")+geom_errorbar(aes(x=6,ymin="+str(round((lower_tau),3))+",ymax="+str(round((upper_tau),3))+",width=.2),color=\"grey\")"
		return UB_s2_string


#==============================================================================
#Main
#==============================================================================


def main(): 

	gene_counts_s1 = defaultdict(list)
	gene_counts_s2 = defaultdict(list)

	for file in args.SE_files:
		readcoverage = str(args.readcoverage)
		passfilt_SE, a_SE, b_SE, c_SE, d_SE = filter_SE(file,readcoverage)

	for file in args.MXE_files:
		readcoverage = str(args.readcoverage)
		passfilt_MXE, a_MXE, b_MXE, c_MXE, d_MXE = filter_MXE(file,readcoverage)



	for i in range(0,len(passfilt_SE)):
		gene = passfilt_SE[i]
		gene_counts_s1[gene].append(sum(a_SE[i]))
		gene_counts_s1[gene].append(sum(b_SE[i]))
		gene_counts_s2[gene].append(sum(c_SE[i]))
		gene_counts_s2[gene].append(sum(d_SE[i]))
		
	for i in range(0,len(passfilt_MXE)):
		gene = passfilt_MXE[i]
		gene_counts_s1[gene].append(sum(a_MXE[i]))
		gene_counts_s1[gene].append(sum(b_MXE[i]))
		gene_counts_s2[gene].append(sum(c_MXE[i]))
		gene_counts_s2[gene].append(sum(d_MXE[i]))

	
	passfilt = passfilt_SE + passfilt_MXE
	a = a_SE + a_MXE
	b = b_SE + b_MXE
	c = c_SE + c_MXE
	d = d_SE + d_MXE

	
	

	a_down, b_down, c_down, d_down = downsample(passfilt, a,b,c,d)

	ijc_s1_norm, sjc_s1_norm, ijc_s2_norm, sjc_s2_norm = normalize(passfilt, a_down, b_down, c_down, d_down)

	#a_dict, b_dict, c_dict, d_dict = make_dictionary(passfilt, ijc_s1_norm, sjc_s1_norm, ijc_s2_norm, sjc_s2_norm)

	gene_max_s1, gene_max_s2 = get_gene_max(passfilt, ijc_s1_norm, sjc_s1_norm, ijc_s2_norm, sjc_s2_norm)

	tau_s1, tau_s2 = get_tau(passfilt, ijc_s1_norm, sjc_s1_norm, ijc_s2_norm, sjc_s2_norm, gene_max_s1, gene_max_s2, gene_counts_s1, gene_counts_s2)


	coordinates = {}
	with open(args.coordinates,"r") as infile:
		for line in infile:
			if not line.startswith("logFC"):
				line = line.rstrip()
				line = line.split("\t")
				gene = line[0]
				bias = line[1]
				coordinates[gene] = bias

	MB_s1 = []
	FB_s1 = []
	UB_s1 = []
	MB_s2 = []
	FB_s2 = []
	UB_s2 = []

	for gene in set(passfilt):
		if gene in coordinates:
			if coordinates[gene] == "MB":
				MB_s1.append(tau_s1[gene])
				MB_s2.append(tau_s2[gene])

			elif coordinates[gene] == "FB":
				FB_s1.append(tau_s1[gene])
				FB_s2.append(tau_s2[gene])

			elif coordinates[gene] == "UB":
				UB_s1.append(tau_s1[gene])
				UB_s2.append(tau_s2[gene])
		# else:
		# 	print "Gene category not found", gene
				#print "UB gene", gene, tau_s1[gene], tau_s2[gene]
		
	print "Number of genes with category info = ", (len(MB_s1)+ len(FB_s1)+ len(UB_s1))

	print "\n"
	print "Total number of genes that pass read depth filter = ", len(set(passfilt))
	print "Total number of junctions that pass read depth filter = ", len(passfilt)

	print "\n"
	MB_s1_string = calculate_stats(MB_s1, "Male-biased expression s1")
	MB_s2_string =calculate_stats(MB_s2, "Male-biased expression s2")
	
	print "\n"
	FB_s1_string =calculate_stats(FB_s1, "Female-biased expression s1")
	FB_s2_string =calculate_stats(FB_s2, "Female-biased expression s2")
	
	print "\n"
	UB_s1_string =calculate_stats(UB_s1, "Unbiased expression s1")
	UB_s2_string =calculate_stats(UB_s2, "Unbiased expression s2")
	

	MB_s1_vals = []
	MB_s2_vals = []
	FB_s1_vals = []
	FB_s2_vals = []
	UB_s1_vals = []
	UB_s2_vals = []

	MB_stat, MB_p = wilcoxon(MB_s1, MB_s2)
	FB_stat, FB_p = wilcoxon(FB_s1, FB_s2)
	UB_stat, UB_p = wilcoxon(UB_s1, UB_s2)

	print "\n"	
	print('MB Statistics=%.3f, p=%.3f' % (MB_stat, MB_p))
	print('FB Statistics=%.3f, p=%.3f' % (FB_stat, FB_p))
	print('UB Statistics=%.3f, p=%.3f' % (UB_stat, UB_p))
	print '\n'
	print('label=c(\"%.3f\",\"%.3f\",\"%.3f\"))' % (MB_p,FB_p,UB_p))
	print '\n'
	print MB_s1_string, MB_s2_string, FB_s1_string, FB_s2_string, UB_s1_string, UB_s2_string
	

if __name__ == "__main__":
	main()

