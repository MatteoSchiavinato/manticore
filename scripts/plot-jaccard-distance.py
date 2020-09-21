#!/usr/bin/env python3

import argparse as ap 
from time import asctime as at
from time import sleep
import pandas as pd 
from plotnine import * 
import sys 
import psutil
import warnings 
import multiprocessing as mp 


warnings.filterwarnings('ignore')


# parser
p = ap.ArgumentParser()
p.add_argument("--combined-file", required=True)
p.add_argument("--gene-annotation")
p.add_argument("--chromosome-lengths", required=True)
p.add_argument("--chromosome")
p.add_argument("--all-chromosomes", action="store_true", default=False)
p.add_argument("--window-size", required=True, type=int)
p.add_argument("--feature", default="whole")
p.add_argument("--basename", default="jacc_analysis")
p.add_argument("--max-cov", type=int, default=50)
p.add_argument("--threads", default=psutil.cpu_count(), type=int)
p.add_argument("--out-dir", default=".")
args = p.parse_args()


if ((args.chromosome==True) and (args.all_chromosomes==True)):
	sys.stderr.write("ERROR: both --chromosome and --all-chromosomes are specified, while only one of them at once should.\n\n")
	sys.exit(1)

# functions 

def filter_raw_table(Raw_jacc_shared_mem, chromosome, feature, max_cov):

	# read in Raw jaccard index table 
	# ...
	# filter keeping only chromosome desired by user 
	# ...
	# filter keeping only feature desired by user 
	# ...
	# keep only relevant statistics (jaccard, mean cov, frac pos) 
	# ...
	# squish coverage to max cov value specified by user 
	# ...
	# normalize everything to 1 for easier plotting (ggplot2 does not allow multi-facet ylim)
	# ...
	# melt table for later plotting 
	# ...
	# return melted table 

	Filt_jacc = Raw_jacc_shared_mem.df[(Raw_jacc['Scaffold'] == chromosome) & (Raw_jacc['Feature'] == feature)]
	Filt_jacc = Filt_jacc.iloc[: , [0,1,5,6,8,9]]
	mask = (Filt_jacc['Mean_cov'] > 50)
	Filt_jacc['Mean_cov'][mask] = 50
	Filt_jacc['W_start'] = Filt_jacc['W_start'] / float(1e6)

	Subgenomes = sorted(list(set(Filt_jacc['Subgenome'].tolist())))

	if len(Subgenomes) > 2:
		sys.stderr.write("There seems to be > 2 subgenomes. This program is made for 2.\n\n")
		sys.exit(1)

	Melt_jacc = pd.melt(	Filt_jacc, id_vars=["W_start", "Subgenome"], 
				value_vars=["Mean_cov", "Frac_pos", "Jaccard"], 
				var_name="Variable", 
				value_name="Value"	)

	Melt_jacc = pd.DataFrame({	'W_start':Melt_jacc['W_start'], 
					'Feature':Melt_jacc['Variable'] + "_" + Melt_jacc['Subgenome'],
					'Value':Melt_jacc['Value']
					})

	keyword = "Jaccard" + "_" + Subgenomes[1]
	Melt_jacc = Melt_jacc[Melt_jacc["Feature"] != keyword]
	Melt_jacc = Melt_jacc.replace("Jaccard" + "_" + Subgenomes[0], "Jaccard").drop_duplicates()
	Melt_jacc.sort_values(by=["W_start", "Feature"], axis='index', inplace=True)

	return (Melt_jacc, Subgenomes, window_size)


def generate_plot(Melt_jacc, Subgenomes, out_dir, basename, chromosome, chrom_length):

	# extract list of subgenomes
 	# ...
	# generate a plot containing five lanes:
	# 1. mean cov SG1
	# 2. mean cov SG2
	# 3. jaccard index
	# 4. frac pos SG1
	# 5. frac pos SG2

	max_w_start = max(Melt_jacc['W_start'].tolist())

	P1 = (
	ggplot(data=Melt_jacc)
	+ theme(plot_title=element_text(family="sans", face="bold", 
					colour='#261A1A', size=11),
	axis_title_x=element_text(	family="sans", face="plain",
					colour='#261A1A', size=8),
	axis_title_y=element_text(	family="sans", face="plain",
					colour='#261A1A', size=8),
	axis_text_x=element_text(family="sans", face="plain", colour='#261A1A', size=8),
	axis_text_y=element_text(family="sans", face="italic", colour='#261A1A', size=8),
	legend_text=element_text(family="sans", face="italic", colour='#000000', size=8),
	legend_title=element_text(family="sans", face="plain", colour='#000000', size=8),
	panel_background=element_rect(fill='#DCDCDC', colour='#DCDCDC'),
	aspect_ratio=0.10)
	+ geom_line(aes(x="W_start", y="Value"))
	+ xlab("Chromosome position [Mbp]")
	+ ylab("")
	+ scale_x_continuous(	breaks=[i for i in range(0, int(round(max_w_start, 0))+1, 10)],
				limits=(0,chrom_length)	)
	+ facet_wrap("~Feature", ncol=1, nrow=5, scales="free")
	)
	
	outfile_base = "{0}/{1}.{2}".format(out_dir, basename, chromosome)
	P1.save(filename="{0}.png".format(outfile_base), format="png", height=4, width=15, dpi=300, units="cm")


def form_windows(Chromosomes, window_size):

	# read scaf lengths file 
	# ...
	# parse each scaffold from 1 to len(scaf)
	# ...
	# save ranges to a dictionary if range longer than window_size
	# ...
	# return Scaf lengths and Windows Dictionary 

	Windows = { scaf : [] for scaf in Chromosomes.keys() if int(Chromosomes[scaf]) >= window_size }
	for scaf in Windows.keys():
		start_pos = 1
		end_pos = int(window_size)
		while (end_pos <= Chromosomes[scaf]):
			Windows[scaf].append( (int(start_pos), int(end_pos)) )
			start_pos += int(window_size)
			end_pos += int(window_size)

	return Windows


def extract_genes(gene_annotation, Windows):

	# read gene annotation into a pandas data frame 
	# ...
	# return genes with their position 

	Genes = pd.read_csv(	gene_annotation, sep="\t", comment="#",
				names=[	"Sequence", "Source", "Feature", "Start", "End", 
					"Score", "Strand", "Phase", "Attributes" ]	).iloc[:,[0,2,3,4]]
	mask = (Genes['Feature'] == "gene")
	Genes = Genes[mask]
	Genes = Genes.iloc[:,[0,2,3]]

	return Genes


def get_gene_densities(scaf, window, Genes_shared_mem, window_size, queue, lock):

	# extract genes within window from data frame 
	# ...
	# count them 
	# ...
	# divide number by window size 
	# ...
	# scale to genes/Mbp
	# ...
	# add value to queue 

	Genes = Genes_shared_mem.df[(Genes_shared_mem.df['Sequence']==str(scaf)) & (((Genes_shared_mem.df['Start']>=int(window[0])) & (Genes_shared_mem.df['Start']<=int(window[1]))) | ((Genes_shared_mem.df['End']>=int(window[0])) & (Genes_shared_mem.df['End']<=int(window[1]))))]

	gene_number = Genes.shape[0]
	conversion_factor = float(1e6) / float(window_size)
	gene_density = float(gene_number) * float(conversion_factor)

	lock.acquire()
	queue.put( (scaf, window[0], gene_density) )
	lock.release()

	return True
	

def run_chromosome_analysis(	Raw_jacc_shared_mem, chromosome, chrom_length, feature, max_cov, 
				out_dir, basename, Densities_shared ):

	# keep only specified feature and chromosome
	# remove irrelevant columns 
	(Melt_jacc, Subgenomes) = filter_raw_table(	Raw_jacc_shared_mem, chromosome, 
							feature, max_cov	)
	# generate plot 
	generate_plot(	Melt_jacc, Subgenomes, out_dir, 
			basename, chromosome, chrom_length	)


def dump_queue(queue):
    
	# Empties all pending items in a queue and returns them in a list
	# ...
	# returns list 
	
	lst = []
	while (queue.qsize() != 0):
		lst.append(queue.get())

	return pd.DataFrame(lst)


# main script 
if __name__ == "__main__":

	sys.stderr.write("\n\n### BEGIN ###\n\n")

	# read input files
	sys.stderr.write("[{0}] Reading input files\n".format(at()))
	Raw_jacc = pd.read_csv(args.combined_file, sep="\t")
	Sequences = list(set(Raw_jacc['Scaffold'].tolist()))
	INPUT = open(args.chromosome_lengths, "r")
	Lsts = [ line.rstrip("\b\r\n").split("\t") for line in INPUT ]
	INPUT.close()
	Chromosomes = { lst[0]:int(lst[1]) for lst in Lsts if lst[0] in Sequences }

	# extract genes 
	sys.stderr.write("[{0}] Extracting genes\n".format(at()))
	Windows = form_windows(Chromosomes, args.window_size)
	Genes = extract_genes(args.gene_annotation, Windows)

	# extract gene densities
	sys.stderr.write("[{0}] Calculating gene densities\n".format(at()))

	manager = mp.Manager()
	Genes_shared_mem = manager.Namespace()
	Genes_shared_mem.df = Genes
	pool = mp.Pool(args.threads)
	lock = manager.Lock()
	queue = manager.Queue()

	Func_args = [ 	(scaf, window, Genes_shared_mem, args.window_size, queue, lock) \
			for scaf in Windows.keys() for window in Windows[scaf] ]

	pool.starmap(func=get_gene_densities, iterable=Func_args)
	pool.close()
	pool.join()

	Gene_densities = dump_queue(queue)
	Gene_densities.columns = ("Sequence", "W_start", "Gene_density")
	Gene_densities.sort_values(by=['Sequence', 'W_start'], axis=0, inplace=True)

	Densities_shared = manager.Namespace()
	Densities_shared.df = Gene_densities

	# generate plot 
	sys.stderr.write("[{0}] Generating picture(s)\n".format(at()))
	Raw_jacc_shared_mem = manager.Namespace()
	Raw_jacc_shared_mem.df = Raw_jacc

	if not args.all_chromosomes:
		Chromosomes = { args.chromosome : int(Chromosomes[args.chromosome]) }

	pool = mp.Pool(args.threads)

	Func_args = [ 	(Raw_jacc_shared_mem, chromosome, round(float(Chromosomes[chromosome])/float(1e6),1), args.feature, args.max_cov, args.out_dir, args.basename, Densities_shared) for chromosome in Chromosomes.keys() ]

	pool.starmap(func=run_chromosome_analysis, iterable=Func_args)

	pool.close()
	pool.join()
	
	sys.stderr.write("\n\n### END ###\n\n")
