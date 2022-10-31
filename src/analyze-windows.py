#!/usr/bin/env python3

### modules ###
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

from time import asctime as at
import numpy as np
import pandas as pd
import dask
import dask.dataframe as dd
from dask.distributed import Client
import dask_distance as ddist
import argparse as ap
import sys
from operator import itemgetter
from itertools import combinations
import shutil

### parser ###
p = ap.ArgumentParser()
p.add_argument("--species-name", required=True)
p.add_argument("--cov-dir", type=str, required=True)
p.add_argument("--cov-names", nargs="*", required=True)
p.add_argument("--output-dir", required=True)
p.add_argument("--beds", nargs="*", default=None)
p.add_argument("--beds-names", nargs="*", default=None)
p.add_argument("--scaf-lengths", required=True)
p.add_argument("--window-size", default=50000, type=int)
p.add_argument("--n-breaks", default=10, type=int)
p.add_argument("--threads", default=1, type=int)
p.add_argument("--min-frac-pos", default=0.01, type=float)
p.add_argument("--min-cov-pos", default=1000, type=int)
p.add_argument("--min-length", default=1000, type=int)
p.add_argument("--min-coverage", default=1, type=int)
p.add_argument("--size-factor", type=float)
args = p.parse_args()

### functions ###

def form_subintervals(scaf_lengths, break_len, window_size):

	# this function takes as input:
	# "scaf lengths": a file containing scaffold names and their length (two columns)
	# "break len": the length of each break, or subinterval, which is args.window_size / args.n_breaks
	# "window size": the value declared with args.window_size, defining the size of each window to operate on

	# read scaf lengths file
	# ...
	# parse each scaffold from 1 to len(scaf)
	# ...
	# save ranges to a dictionary if range longer than window_size
	# ...
	# return Scaf lengths and All_subintervals Dictionary

	INPUT = open(scaf_lengths, "r")
	Lines = [ line.rstrip("\b\r\n").split("\t") for line in INPUT ]
	Scaf_lengths = { str(lst[0]) : int(lst[1]) for lst in Lines }
	INPUT.close()

	All_subintervals = { scaf : [] for scaf in Scaf_lengths.keys() if int(Scaf_lengths[scaf]) >= window_size }
	for scaf in All_subintervals.keys():
		scaf_len = Scaf_lengths[scaf]
		start_pos = 1
		end_pos = int(break_len)
		# the end of a scaffold is not used
		# because the while loop breaks when end pos > scaf_length
		# Future: try to keep that too
		while (end_pos <= Scaf_lengths[scaf]):
			All_subintervals[scaf].append( (int(start_pos), int(end_pos)) )
			start_pos += int(break_len)
			end_pos += int(break_len)

	return All_subintervals


def compute_real_lengths(All_subintervals, bed_file, window_size):

	# this function reads the intervals provided in a BED file
	# these intervals define the regions of the scaffolds where the information must be analysed
	# the positions outside of these regions are not used for the analysis
	# Note: if no bed file is provided by the user, these regions are not computed and the whole sequence is used

	# The regions are read line by line from the BED file
	# Each region is then compared to the Subintervals computed in the previous function
	# The overlapping positions between each break (or subinterval) and each annotated feature is termed as "real length"
	# Real length indicates the "real length" of the break, i.e. the length that was used to perform the analysis
	# Real length must be <= break length
	# If a feature overlaps a break but extends beyond it:
	# the overhang is added to the next break (see last part of the function)
	# Note: this overhang addition is performed until there is an overhang
	# this means that if the overhang spans the next 3 breaks, all of them will receive a contribution from the overhang
	# (and not just the next break)

	# initiate real lengths dictionary
	# ...
	# for line in bed file
	# ...
	# find suitable subinterval (if existing) -> match with All_subintervals dictionary
	# ...
	# compute overlap between feature and subinterval -> real length
	# ...
	# if overhang is present: add it to next subinterval until there are valid subintervals
	# ...
	# return real lengths dictionary and number of skipped lines

	Real_lengths = { scaf : { pos_range : 0 for pos_range in All_subintervals[scaf] } \
			for scaf in All_subintervals }

	INPUT = open(bed_file, "r")
	k=0
	old_scaf = ""

	for line in INPUT:
		lst = line.rstrip("\n\r\b").split("\t")
		(scaf, start, end) = (lst[0], int(lst[1])+1, int(lst[2]))
		if scaf in All_subintervals.keys():

			# check last subinterval: is start bigger than that? then don't bother
			if not (start > All_subintervals[scaf][-1][1]):
				# restart the counter if scaffold changed
				if scaf != old_scaf:
					k=0
				# add +1 to counter if start is larger than end of current subinterval
				while (	(start > All_subintervals[scaf][k][1]) and \
					(k < len(All_subintervals[scaf])-1)):
					k+=1

				# compute overlap between features and subintervals
				# sum the overlap to the real length of the subinterval
				overlap_start = max(start, All_subintervals[scaf][k][0])
				overlap_end = min(end, All_subintervals[scaf][k][1])
				overlap_length = overlap_end - overlap_start + 1
				pos_range = All_subintervals[scaf][k]
				Real_lengths[scaf][pos_range] += overlap_length

				# if after summing up there is still some "feature" that overhangs...
				if (end > All_subintervals[scaf][k][1]):
					j=1
					overhang_length = end - All_subintervals[scaf][k][1]

					# while there are subintervals after the current one
					while ((k+j <= len(All_subintervals[scaf])-1) and (overhang_length > 0)):
						next_pos_range = All_subintervals[scaf][k+j]
						Real_lengths[scaf][next_pos_range] += min(overhang_length, window_size)
						j+=1
						overhang_length = int(overhang_length) - int(window_size)
						# if this was the last iteration, the number drops < 0
						# I set it at 0 to stop the while loop
						if overhang_length < 0:
							overhang_length = 0

		old_scaf = scaf

	INPUT.close()

	return Real_lengths


def generate_dataframes(args, region):

	####################
	# define input files
	####################

	# input files are declared here
	# the user must know that this script is not a standalone script so it should not be run independently
	# however, if you change these two files with the ones that you want to use, it will work
	infile_1 = "{0}/{1}.{2}.{3}x.depth".format(args.cov_dir, args.cov_names[0], region, args.min_coverage)
	infile_2 = "{0}/{1}.{2}.{3}x.depth".format(args.cov_dir, args.cov_names[1], region, args.min_coverage)

	#####################
	# generate data frame
	#####################

	# the two coverage files (one per parent) are read as DASK dataframes
	# dask has the same syntax as pandas but lazy-loads the dataframes
	# this means that, unless a computation is triggered, operations to be performed on the dataframe are accumulated
	# once the compute() method is called, the computations are executed all at once
	# they are executed making use of the multi-threaded potential of dask
	df1 = dd.read_csv(infile_1, sep="\t", names=["Sequence", "Position", "Coverage"], dtype={"Sequence":"category", "Position":"uint32", "Coverage":"uint32"})
	df1 = df1.repartition(npartitions=args.threads)
	df2 = dd.read_csv(infile_2, sep="\t", names=["Sequence", "Position", "Coverage"], dtype={"Sequence":"category", "Position":"uint32", "Coverage":"uint32"})
	df2 = df2.repartition(npartitions=args.threads)

	return (df1, df2)


def filter_dataframes(df1, df2, args, All_subintervals, Scaf_lengths):

	# these variables are initialized
	# a list of all the scaffolds present in the coverage files and in the All_subintervals dictionary is generated
	# the list is an intersection

	# seqnum is the counter of how many sequences have been processed
	# totseq is the total number of basepairs of the selected scaffolds
	# the counter will advance in a way that is proportional to the number of basepairs computed
	# Example:
	# If the genome has two chromosomes, one with 8 windows and one with 2 windows
	# the first chromosome will advance the % counter by 80%, and the second one by 20%

	# I count them from "All_subintervals" because the keys of that dictionary are the only sequences that will be processed
	# This is because, when I create the "All_subintervals" dictionary, I make a selection of sequences that are at least as long as the window size
	# It would be pointless to process the other sequences since they wouldn't be able to return a single window
	# I then iterate the main code operations for each sequence of those included in "All_subintervals"
	# In assemblies that are highly fragmented, this can result in a very drastic reduction of computational times
	Sub_scafs = list(All_subintervals.keys())
	Cov_scafs_1 = list(set(df1["Sequence"].compute(num_workers=args.threads).tolist()))
	Cov_scafs_2 = list(set(df2["Sequence"].compute(num_workers=args.threads).tolist()))
	Scafs = list(set(Sub_scafs) & (set(Cov_scafs_1) | set(Cov_scafs_2)))

	Scaf_lengths_sub = { scaf : Scaf_lengths[scaf] for scaf in Scafs }

	seqnum = 0
	totseq = sum(Scaf_lengths_sub.values())
	error_message = "[{0}] {1} sequences could produce windows of {2} bp size\n".format(at(), len(Scafs), args.window_size)

	# a reduction of the size of "df1" and "df2" is made, based on the scffolds that are in All_subintervals
	# this means that only scaffolds containing at least args.window_size base pairs are retained
	df1 = df1[df1["Sequence"].isin(Scafs)].repartition(npartitions=args.threads)
	df1["Sequence"] = df1["Sequence"].cat.set_categories(Scafs)
	df2 = df2[df2["Sequence"].isin(Scafs)].repartition(npartitions=args.threads)
	df2["Sequence"] = df2["Sequence"].cat.set_categories(Scafs)

	return (df1, df2, error_message)


def estimate_size_factor(args, df1, df2):

	# a user can define a size factor between the sequencing libraries used to produce the coverage files
	# this size factor is the ratio between the expected mean coverages produced by the libraries
	# if this is not declared, the script computes it for you
	# the calculation is done using the observed coverage
	# the coverage column of a coverage file is extracted and the average is computed
	# then the MEAN coverage in parent 1 is divided by the one in parent 2
	# this returns the size factor

	# Note:
	# this factor is used later -> the coverage in df2 (parent2) is multiplied by the size factor
	# This makes it comparable with parent 1

	# Note #2:
	# the size factor is computed with the following assumptions:
	# - reads from a parent will mostly map to regions of the genome belonging to that parent
	# - the produced mean coverage X will be an approximation of the expected coverage
	# - due to evolution, we expected obs(X) < exp(X) (mutations and indels will hinder mapping of some reads, reducing X)

	# Note #3:
	# Biased fractionation should not affect this value too much
	# In fact, this value DOES NOT depend on the genome (or subgenome) size
	# It is only computed using the positions covered by each subgenome
	# That means that, even if one subgenome is N times larger than the other, the mean coverage does not depend on this

	if not args.size_factor:

		lib_cov_exp_1 = df1["Coverage"].mean().compute(num_workers=args.threads)
		lib_cov_exp_2 = df2["Coverage"].mean().compute(num_workers=args.threads)
		sizeFactor = float(lib_cov_exp_1) / float(lib_cov_exp_2)

		error_message = "[{0}] A size factor of {1} has been estimated between {2} and {3}\n".format(at(), sizeFactor, args.cov_names[0], args.cov_names[1])

	else:

		sizeFactor = args.size_factor
		error_message = "[{0}] A size factor of {1} has been declared between {2} and {3}\n".format(at(), sizeFactor, args.cov_names[0], args.cov_names[1])

	return (sizeFactor, error_message)


def add_subinterval(df1, df2, All_subintervals, Real_lengths, break_len, region, args):

	# the subinterval and the real length columns are added
	# these columns are needed for the calculation of the subinterval metrics (see later)
	# they are added using simple lambda functions
	def add_subinterval_and_remove_end_positions(df1, df2, break_len, All_subintervals):

		# SUBINTERVAL
		# this is computed using a number conversion
		# x is the position in the scaffold, in 1-based coordinates
		# (x-1) is the position in 0-based coordinates

		# (x-1) % break_len returns the break in 0-based coordinates
		# (x-1) % break_len +1 returns the break in 1-based coordinates, which is what we ultimately need
		# example:
		# position is x=5751 , break_len = 5000
		# (x-1) % break_len +1 = 2
		# In fact, we are in the second subinterval (which goes from 5001 to 10000)
		df1['Subinterval'] = df1['Position'].apply(lambda x : int((x-1) - ((x-1) % break_len) + 1), meta="uint32")
		df1['Subinterval'] = df1['Subinterval'].astype("uint32")
		df2['Subinterval'] = df2['Position'].apply(lambda x : int((x-1) - ((x-1) % break_len) + 1), meta="uint32")
		df2['Subinterval'] = df2['Subinterval'].astype("uint32")

		# removal of positions that are outside subintervals
		# the end of sequences cannot make an entire window, so they are chomped
		# First, the subinterval column is added
		# Second, the Real length column is initialized
		# Third, remove all those positions that go beyond the scaffold limit
		# Fourth, add the real length of a certain subinterval
		df1 = df1[df1['Position'] <= df1['Sequence'].apply(lambda x : All_subintervals[x][-1][1], meta="uint32").astype("uint32")]
		df2 = df2[df2['Position'] <= df2['Sequence'].apply(lambda x : All_subintervals[x][-1][1], meta="uint32").astype("uint32")]

		# reorder columns
		# Sequence | Subinterval | Coverage
		df1 = df1.loc[:,["Sequence", "Subinterval", "Coverage"]]
		df2 = df2.loc[:,["Sequence", "Subinterval", "Coverage"]]

		return (df1, df2)

	(df1, df2) = add_subinterval_and_remove_end_positions(df1, df2, break_len, All_subintervals)

	return (df1, df2)


def process_subintervals(df1, df2, Real_lengths, region, break_len):

	# metrics are computed for each subinterval (see function for details)
	# the pandas (or dask) groupby() method "groups" lines of a dataframe by a certain category that is defined in one of the columns
	# in this case, we previously declared the subinterval column
	# each different subinterval in that column will represent a group name
	# all the lines carrying that subinterval annotation will be its group
	# all the metrics are computed on that group and returned
	# the outcome is a dataframe where line names are subintervals and each column is a metric

	# nested:
	def compute_metrics(x, Real_lengths, region, break_len):

		# The metrics for each subinterval are computed
		# These metrics are not the final metrics, but are a necessary mid-step
		# In fact, the metrics for each window will be computed starting from the metrics of the window's subintervals

		# REAL LENGTH
		# In this step the "Real_length" of the subinterval is computed as the mean of the real lengths annotated in the Real_length column
		# Although it might look counter-intuitive, since a subinterval has only ONE real length, the reason is:
		# We are still at a stage where the dataframe has the Sequence, Position and Coverage information for each POSITION
		# To perform a pandas groupby() operation, the program has added the "Subinterval" column and the corresponding Real length
		# This allows us to use the pd.groupby() function, grouping by subinterval
		# As a trade-off for that, I am computing the real length as a mean of the real length column
		# but the real length column contains the same value for each line, that is, the real length of the subinterval
		# In other words, this is not really a meaningful mean, because is the mean of a series of identical values

		# COV POS
		# This metric represents the number of covered position within a break
		# It is computed just by counting the number of positions annotated for each break

		# FRAC POS
		# This metric represents the fraction of a break that is covered
		# It is computed by dividing the Cov pos by the Real length

		# MEAN COV
		# This is the final metric and is computed by averaging the values in the Coverage column

		# compute real length
		# ...
		# sum total covered positions
		# ...
		# divide covered positions by real length
		# ...
		# average coverage per position
		# ...
		# return a Series that represents a subinterval

		d = {}
		scaf = x["Sequence"].iloc[0]
		start = x["Subinterval"].iloc[0]
		end = start + break_len -1
		subinterval = (start, end)
		real_length = Real_lengths[region][scaf][subinterval]
		d["Cov_pos"] = x.shape[0]
		try:
			d["Frac_pos"] = float(d["Cov_pos"]) / float(real_length)
		except ZeroDivisionError:
			d["Frac_pos"] = float(0)

		d["Mean_cov"] = x["Coverage"].mean()
		return pd.Series(d, index=["Cov_pos", "Frac_pos", "Mean_cov"])

	df1 = df1.groupby(["Sequence", "Subinterval"]).apply(compute_metrics, meta={"Cov_pos":"uint32", "Frac_pos":"float32", "Mean_cov":"float32"}, Real_lengths=Real_lengths, region=region, break_len=break_len)
	df1 = df1.dropna(how="all")
	df2 = df2.groupby(["Sequence", "Subinterval"]).apply(compute_metrics, meta={"Cov_pos":"uint32", "Frac_pos":"float32", "Mean_cov":"float32"}, Real_lengths=Real_lengths, region=region, break_len=break_len)
	df2 = df2.dropna(how="all")
	df1.columns = ["Cov_pos_1", "Frac_pos_1", "Mean_cov_1"]
	df2.columns = ["Cov_pos_2", "Frac_pos_2", "Mean_cov_2"]

	return (df1, df2)


def merge_dataframes(df1, df2, Real_lengths, region, break_len, args):

	####################
	# compute dataframes
	# merge dataframes
	# bring them back to dask
	#########################

	# computation of both dataframes is triggered via .compute()
	# the larger the scaffold, the more parallel operations will be performed
	# this because each thread is a subinterval
	# the computation returns a pandas dataframe for each parent
	# the two parental pd dataframes are concatenated by column
	# the index is reset
	# the final concatenated df is passed again to dask, to avoid keeping it in memory for no reason
	# substitute NaN with zeros
	# then the types of each column are re-declared for later analyses
	df1 = df1.compute(num_workers=args.threads).sort_index()
	df2 = df2.compute(num_workers=args.threads).sort_index()
	df = df1.join(df2, on=["Sequence", "Subinterval"], how="outer").sort_index()
	df = df.fillna(value=0)
	df = df.reset_index()
	df = dd.from_pandas(df, npartitions=args.threads)

	# add real length column
	def get_real_length(x, Real_lengths, region, break_len):

		scaf = x["Sequence"]
		start = x["Subinterval"]
		end = start + break_len -1
		subinterval = (start, end)
		real_length = int(Real_lengths[region][scaf][subinterval])
		return real_length

	df["Real_length"] = df[["Sequence", "Subinterval"]].apply(get_real_length, axis=1, meta="uint32", Real_lengths=Real_lengths, region=region, break_len=break_len)
	df = df.loc[:,["Sequence", "Subinterval", "Real_length", "Cov_pos_1", "Frac_pos_1", "Mean_cov_1", "Cov_pos_2", "Frac_pos_2", "Mean_cov_2"]]
	df = df.astype({"Subinterval":"uint32", "Real_length":"uint32", "Cov_pos_1":"uint32", "Frac_pos_1":"float32", "Mean_cov_1":"float32", "Cov_pos_2":"uint32", "Frac_pos_2":"float32", "Mean_cov_2":"float32"})

	return df


def analyze_windows(df, args, region, break_len):

	####################################################
	# generate jaccard index arrays to be computed later
	# this function contains TWO HARDCODED values
	# these might become options in the Future
	# for now, they are like This
	# the user doesn't know any better, anyway, on how to set them
	def generate_jaccard_arrays(df):

		###########################
		# add jaccard array columns
		# they contain a value of True / False depending on if the subinterval is covered by a parent
		# this value is the modified as follows:
		# if other mean cov is at more than 25% distance (which means ratio is < 0.75)
		# and if frac pos is more than 10%
		# turn to True
		# --------------------------
		# these values are HARDCODED
		# the user might want to change them but they are "too much" in the options
		###########################################################################

		MIN_BREAK_MEAN_COV_RATIO = 0.75
		MIN_BREAK_FRAC_POS = 0.1

		df['Jacc_array_1'] = False
		df['Jacc_array_2'] = False

		# turn to true those subintervals which have the sufficient fraction of covered positions
		condition = (df["Frac_pos_1"] >= MIN_BREAK_FRAC_POS)
		df["Jacc_array_1"] = df["Jacc_array_1"].mask(condition, True)
		condition = (df["Frac_pos_2"] >= MIN_BREAK_FRAC_POS)
		df["Jacc_array_2"] = df["Jacc_array_2"].mask(condition, True)

		# revert to false those whose coverage falls below 75% of the other
		condition = (df["Mean_cov_1"] / df["Mean_cov_2"] < MIN_BREAK_MEAN_COV_RATIO)
		df["Jacc_array_1"] = df["Jacc_array_1"].mask(condition, False)
		condition = (df["Mean_cov_2"] / df["Mean_cov_1"] < MIN_BREAK_MEAN_COV_RATIO)
		df["Jacc_array_2"] = df["Jacc_array_2"].mask(condition, False)

		return df

	# add window column for later grouping
	# the window line is added similarly to how the subinterval line was added
	# then, column order is rearranged
	df = generate_jaccard_arrays(df)
	df['Window'] = df["Subinterval"].apply(lambda x : x-(x % args.window_size)+1, meta="uint32")
	df = df[["Sequence", "Real_length", "Window", "Cov_pos_1", "Frac_pos_1", "Mean_cov_1", "Cov_pos_2", "Frac_pos_2", "Mean_cov_2", "Jacc_array_1", "Jacc_array_2"]]

	################################
	# calculate window-based metrics

	# First, dataframe is grouped by sequence
	# Second, dataframe is grouped by windows
	# Third, the function is invoked
	# metrics for each window are computed (see function for details)
	# In this step, the .compute() method is invoked
	# From here on, I use pandas and not dask because the dataframe is small enough to be handled with pandas directly
	# the index is then re-set
	# a series of columns are added / filled-in:
	# scaffold, feature, subgenome (just the column, not the content)
	# column order is rearranged
	def compute_window_metrics(x):

		# This function is an evolution of the previous one
		# It computes the same metrics but it does so for parent 1 and 2 (which are now combined in a single dataframe)
		# Note:
		# this time the Real length is computed as a sum() and not as a mean(). That is because I am summing up the real length of each subinterval
		# the same is done for Cov pos: I don't count the number of lines but I sum up the values in the cov pos field
		# these values were previously computed for each subinterval

		# COV POS, FRAC POS, MEAN COV
		# These metrics represent the same as they do in the Subintervals
		# The only difference is the "mean_cov"
		# IN this case, I am computing it by taking only the subintervals that have a coverage
		# This way, I ensure that the mean cov of the window is not below the specified minimum coverage passed from command line
		# If it is below, it is 0

		# v1 and v2 ARRAYS
		# I create two numpy arrays composed of boolean [True, False] values
		# These arrays are composed by reading the "Jacc_array_*" column of parents 1 and 2
		# These columns are generated outside of the function
		# They contain a True value if the parent has coverage in that specific subinterval
		# They contain a False value if it doesn't
		# Having coverage is not sufficient: there are some thresholds set to assign a True value (see main code block)

		# UNION
		# The script performs a logical union of the v1 and v2 arrays
		# The union includes all subintervals whereby there is a True in either one of the parents (logical OR)

		# INTERSECTION
		# This works like the union but the intersection is made
		# The intersection includes all subintervals whereby there is a True for BOTH parents (logical AND)

		# UNIQUE
		# This is the number of subintervals that are covered by one parent and NOT by the other
		# It is performed using standard python logical operators
		# ~v2 returns the opposite of v2
		# v1 & ~v2 returns v1 AND the opposite of v2, that is, where only v1 is True
		# same for v2 later

		# JACCARD
		# The Jaccard distance is computed dividing the intersection by the union
		# The closer the intersection gets to the union, the closer the J gets to 1
		# A J=1 means complete overlap between datasets
		# A J=0 means no overlap between data sets
		# In biological terms:
		# J=1 means that there is contribution from both parents in each subinterval
		# J=0 means that, in each subinterval, there is contribution from only one of the two parents
		# Note:
		# The size of the subintervals determins how much "intermixing" you observe
		# In fact, smaller breaks will likely lower Jaccard indexes
		# because it is just less likely to have contribution from both parents in a smaller interval

		# sum real lengths of subintervals
		# ...
		# sum covered positions of subintervals
		# ...
		# average fractions of covered positions
		# ...
		# average mean coverages
		# ...
		# calculate a jaccard distance
		# ...
		# return a Series that represents a window

		d = {}
		d["Real_length"] = x["Real_length"].sum()
		# parent 1
		d["Cov_pos_1"] = x["Cov_pos_1"].sum()
		d["Frac_pos_1"] = x["Frac_pos_1"].mean()
		# assign zero if there are no subintervals
		if x[x["Mean_cov_1"] > 0].shape[0] > 0:
			d["Mean_cov_1"] = x[x["Mean_cov_1"] > 0]["Mean_cov_1"].mean()
		else:
			d["Mean_cov_1"] = 0
		# parent 2
		d["Cov_pos_2"] = x["Cov_pos_2"].sum()
		d["Frac_pos_2"] = x["Frac_pos_2"].mean()
		# assign zero if there are no subintervals
		if x[x["Mean_cov_2"] > 0].shape[0] > 0:
			d["Mean_cov_2"] = x[x["Mean_cov_2"] > 0]["Mean_cov_2"].mean()
		else:
			d["Mean_cov_2"] = 0

		v1 = np.array(x["Jacc_array_1"])
		v2 = np.array(x["Jacc_array_2"])
		d["Union"] = np.logical_or(v1, v2).sum()
		d["Intersection"] = np.logical_and(v1, v2).sum()
		d["Unique_1"] = np.array(v1 & ~v2).sum()
		d["Unique_2"] = np.array(v2 & ~v1).sum()
		if (d["Union"] > 0):
			d["Jaccard"] = float(d["Intersection"]) / float(d["Union"])
		else:
			d["Jaccard"] = -0.1

		return pd.Series(d, index=[	"Real_length",
									"Cov_pos_1", "Frac_pos_1", "Mean_cov_1",
									"Cov_pos_2", "Frac_pos_2", "Mean_cov_2",
									"Jaccard", "Union", "Intersection", "Unique_1", "Unique_2"])

	# execute window analysis
	df = df.groupby(["Sequence", "Window"]).apply(	compute_window_metrics,
													meta={	"Real_length":"uint32",
															"Cov_pos_1":"uint32", "Frac_pos_1":"float32", "Mean_cov_1":"float32",
															"Cov_pos_2":"uint32", "Frac_pos_2":"float32", "Mean_cov_2":"float32",
															"Jaccard":"float32", "Union":"int64", "Intersection":"int64", "Unique_1":"int64", "Unique_2":"int64"})

	df = df.dropna(how="all")
	df = df.reset_index()
	df = df.compute(num_workers=args.threads).sort_values(by=["Sequence", "Window"])

	# change Window to W_start and W_end
	df['W_start'] = df['Window']
	df['W_end'] = df['W_start'] + args.window_size -1
	df = df.loc[:,[	"Sequence", "W_start", "W_end", "Real_length",
				 	"Cov_pos_1", "Frac_pos_1", "Mean_cov_1",
					"Cov_pos_2", "Frac_pos_2", "Mean_cov_2",
					"Jaccard", "Union", "Intersection", "Unique_1", "Unique_2"]]

	df = df.sort_values(by=["Sequence", "W_start"])
	return df


def complete_data_frame(df, region, break_len, args):

	df["Feature"] = region
	df["Subgenome"] = "NA"
	df = df[["Sequence", "W_start", "W_end", "Feature", "Real_length", "Jaccard", "Subgenome", "Cov_pos_1", "Frac_pos_1", "Mean_cov_1", "Cov_pos_2", "Frac_pos_2", "Mean_cov_2", "Union", "Intersection", "Unique_1", "Unique_2"]]

	# Intersection, Union, Unique_1/2 and Uncovered values are adapted
	# At this point, these columns have the number of intersection, union and uncovered BREAKS
	# What we are interested in, however, is the number of POSITIONS
	# While we could compute the position values using the raw coverage file, we did the whole analysis using breaks as units
	# that means that it is actually more accurate to keep doing so
	# This means that, to compute Union, Intersection and Uncovered I just use the break_len and the window_size
	# Note:
	# In case the entire genome is covered this might return total values that are slightly larger than the assembly size
	# This is because, in reality, end-of-scaffold short windows are not used since they are too short to fit the --window-size
	# Later in the main Manticore script, this is handled by taking the:
	# "min(scaffold_size, <intersection | union | uncovered>)"
	df["Intersection"] = df["Intersection"] * break_len
	df["Union"] = df["Union"] * break_len
	df["Unique_1"] = df["Unique_1"] * break_len
	df["Unique_2"] = df["Unique_2"] * break_len
	df["Uncovered"] = (df["Union"] - args.window_size) * -1

	# set data types
	df = df.astype({"Sequence":"object", "W_start":"uint32", "W_end":"uint32", "Feature":"object", "Real_length":"uint32",
					"Jaccard":"float32", "Subgenome":"object",
					"Cov_pos_1":"uint32", "Frac_pos_1":"float32", "Mean_cov_1":"float32",
					"Cov_pos_2":"uint32", "Frac_pos_2":"float32", "Mean_cov_2":"float32",
					"Union":"uint32", "Intersection":"uint32", "Unique_1":"uint32", "Unique_2":"uint32", "Uncovered":"uint32"})

	# set [ J = -0.1 ] for lines not fitting requirements
	# if sum of fraction of covered positions is below --min-frac-pos
	# or if sum of covered positions is below --min-cov-pos
	# or if real length below --min-length
	# reject the computed jaccard index
	# by applying an out-of-bound value (-0.1)
	mask = (df[["Frac_pos_1", "Frac_pos_2"]].sum(axis=1) < float(args.min_frac_pos)) | (df[["Cov_pos_1", "Cov_pos_2"]].sum(axis=1) < float(args.min_cov_pos)) | (df["Real_length"] < float(args.min_length))
	df[mask]["Jaccard"] = -0.1

	# df1 and df2 are generated again by splitting the main df
	df1 = df.iloc[:,[0,1,2,3,4,5,6,7,8,9,13,14,15,17]]
	df2 = df.iloc[:,[0,1,2,3,4,5,6,10,11,12,13,14,16,17]]

	# rename columns
	Final_columns = ["Sequence", "W_start", "W_end", "Feature", "Real_length", "Jaccard", "Subgenome", "Cov_pos", "Frac_pos", "Mean_cov", "Union", "Intersection", "Unique", "Uncovered"]
	df1.columns = Final_columns
	df2.columns = Final_columns

	# recast data types
	df1 = df1.astype({"Sequence":"object", "W_start":"uint32", "W_end":"uint32", "Feature":"object", "Real_length":"uint32", "Jaccard":"float32", "Subgenome":"object", "Cov_pos":"uint32", "Frac_pos":"float32", "Mean_cov":"float32", "Union":"uint32", "Intersection":"uint32", "Unique":"uint32", "Uncovered":"uint32"})
	df2 = df2.astype({"Sequence":"object", "W_start":"uint32", "W_end":"uint32", "Feature":"object", "Real_length":"uint32", "Jaccard":"float32", "Subgenome":"object", "Cov_pos":"uint32", "Frac_pos":"float32", "Mean_cov":"float32", "Union":"uint32", "Intersection":"uint32", "Unique":"uint32", "Uncovered":"uint32"})

	# assign subgenome
	df1["Subgenome"] = args.cov_names[0]
	df2["Subgenome"] = args.cov_names[1]

	# concatenate
	# if the final df does not exist, initiate it with the first scaffold
	# otherwise just append the new scaffold df to the main, final df
	df = pd.concat([df1, df2], axis=0)

	return df


### MAIN PROGRAM ###

if __name__ == "__main__":
	sys.stderr.write("\n[{0}] ### BEGIN ###\n\n".format(at()))

	#################
	# session startup
	#################

	# the command run by the user is printed to STDERR
	# (in case the user has to re-run it quickly or forgot which parameters were used)
	command = " ".join(sys.argv)
	sys.stderr.write("{0}\n\n".format(command))

	# the break length is computed
	# if the computation returns a float (that is, if there is a %) then this is not right
	# subintervals (i.e. breaks) must be integers
	# so if there is a float, the script dies
	if int(args.window_size) % int(args.n_breaks) != 0:
		sys.exit("\nERROR: --window-size divided by --n-breaks needs to return an integer\n\n")
	else:
		break_len = int(float(args.window_size) / float(args.n_breaks))

	###################
	# form subintervals
	# read scaf lengths
	###################

	sys.stderr.write("[{0}] Calculating windows and position ranges ... ".format(at()))

	# subintervals are computed (see function for details)
	All_subintervals = form_subintervals( args.scaf_lengths, break_len, args.window_size )
	INPUT = open(args.scaf_lengths, "r")
	Lsts = [ line.rstrip("\b\r\n").split("\t") for line in INPUT ]
	INPUT.close()
	Scaf_lengths = { lst[0]:int(lst[1]) for lst in Lsts }

	sys.stderr.write("DONE\n")

	######################
	# compute real lengths
	# using provided bed files
	# and provided region names
	# real length is computed for the BREAK and not for the entire WINDOW
	# WINDOW is used for Jaccard index calculation
	##############################################

	sys.stderr.write("[{0}] Computing length of features within windows ... ".format(at()))

	# real lengths are computed (See function for details)
	if (args.beds == None) and (args.beds_names == None):
		Names = ["whole"]
		Real_lengths = { "whole" : { scaf : { subinterval : break_len \
				for subinterval in All_subintervals[scaf] } \
				for scaf in All_subintervals.keys() } }
	else:
		Names = [i for i in args.beds_names]
		Beds = { args.beds_names[i] : args.beds[i] for i in range(0, len(args.beds)) }
		Real_lengths = {}
		for region in Beds.keys():
			Region_real_lengths = compute_real_lengths(	All_subintervals,
									Beds[region],
									break_len )
			Real_lengths[region] = Region_real_lengths

	sys.stderr.write("DONE\n")

	###################
	# coverage analysis
	###################

	# this will become a list of dataframes
	# each dataframe corresponds to the coverage data in a certain region (i.e. feature, like CDS or nonrep)
	Dataframes = []

	# the analysis is performed region by region
	# this means that if --args.beds was declared, the script proceeds one bed file at a time
	# for each bed file, parents are analysed simultaneously

	for region in Names:
	# for region in args.beds_names:

		sys.stderr.write("[{0}] Working on region: {1}\n".format(at(), region))

		######################
		# generate data frames

		sys.stderr.write("[{0}] Generating data frames ... ".format(at()))
		(df1, df2) = generate_dataframes(args, region)
		sys.stderr.write("DONE\n")

		######################
		# estimate size factor

		sys.stderr.write("[{0}] Estimating size factor between coverage profiles ... ".format(at()))
		(sizeFactor, error_message) = estimate_size_factor(args, df1, df2)
		sys.stderr.write("DONE\n")
		sys.stderr.write(error_message)

		###############################
		# keep only necessary sequences

		sys.stderr.write("[{0}] Filtering data frames ... ".format(at()))
		(df1, df2, error_message) = filter_dataframes(df1, df2, args, All_subintervals, Scaf_lengths)
		sys.stderr.write("DONE\n")
		sys.stderr.write(error_message)

		########################
		# add subinterval column

		sys.stderr.write("[{0}] Extracting valid positions ... ".format(at()))
		(df1, df2) = add_subinterval(df1, df2, All_subintervals, Real_lengths, break_len, region, args)
		sys.stderr.write("DONE\n")

		#############################
		# compute subinterval metrics

		sys.stderr.write("[{0}] Computing metrics at a break level ... ".format(at()))
		(df1, df2) = process_subintervals(df1, df2, Real_lengths, region, break_len)
		sys.stderr.write("DONE\n")

		###################
		# apply size factor

		sys.stderr.write("[{0}] Applying size factor to parent 2 ... ".format(at()))
		# the size factor is applied (see previous description)
		# it is just used in a simple multiplication
		df2['Mean_cov_2'] = df2['Mean_cov_2'].apply(lambda x : x * sizeFactor, meta="float32")
		df1 = df1.astype({"Cov_pos_1":"uint32", "Frac_pos_1":"float32", "Mean_cov_1":"float32"})
		df2 = df2.astype({"Cov_pos_2":"uint32", "Frac_pos_2":"float32", "Mean_cov_2":"float32"})
		sys.stderr.write("DONE\n")

		###################################
		# merge the two dataframes into one

		sys.stderr.write("[{0}] Integrating info from both parents in one data frame ... ".format(at()))
		df = merge_dataframes(df1, df2, Real_lengths, region, break_len, args)
		sys.stderr.write("DONE\n")

		########################
		# compute window metrics

		sys.stderr.write("[{0}] Computing window-based metrics ... ".format(at()))
		df = analyze_windows(df, args, region, break_len)
		sys.stderr.write("DONE\n")

		#####################
		# finalize data frame
		sys.stderr.write("[{0}] Finalizing dataframe ... ".format(at()))
		df = complete_data_frame(df, region, break_len, args)
		sys.stderr.write("DONE\n")

		#####################
		# append to main list

		Dataframes.append(df)


	##############################
	# concatenating all dataframes
	##############################

	main_df = pd.concat(Dataframes)

	###################
	# writing to output
	###################

	# adjust values
	# frac pos is brought from [0,1] to [0,100] for human readability in %
	# float values are rounded to the 2nd digit, again for human readability
	main_df["Frac_pos"] = main_df["Frac_pos"] * float(100)
	main_df = main_df.round({"Jaccard":2, "Frac_pos":2, "Mean_cov":2})

	# sort
	main_df.sort_values(by=["Sequence", "W_start", "W_end", "Feature", "Subgenome"], inplace=True)

	# main file
	main_df.to_csv("{0}/{1}.combined.results.txt".format(args.output_dir, args.species_name), sep="\t", index=False)

	# genome-specific files
	for region in args.beds_names:
		for name in args.cov_names:
			outfile = "{0}/{1}.{2}.{3}.txt".format(args.output_dir, args.species_name, name, region)
			df_sub = main_df[(main_df["Subgenome"]==name) & (main_df["Feature"]==region)]
			df_sub.to_csv(outfile, sep="\t", index=False)

	sys.stderr.write("\n[{0}] ### END ###\n\n".format(at()))
