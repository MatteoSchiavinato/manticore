#!/usr/bin/env python3

import argparse as ap
import sys
import os 
import os.path 
import multiprocessing as mp 
import time
from time import asctime as at
from operator import itemgetter
from shutil import copyfile 
import subprocess

start_time = time.time()

VERSION = '''

--------------------------------------------------------
Matteo Schiavinato, Alexandrina Bodrug
BOKU - University of Natural Resources and Life Sciences
Vienna (AT) 
--------------------------------------------------------

'''

if len(sys.argv) < 2:
	sys.argv.append("-h")

if sys.argv[1] in ["-h", "-help", "--help", "getopt", "usage"]:

	sys.exit('''

{0}

Sometimes you want to re-run our pipeline at a different minimum coverage, higher than the one 
you originally used, or simply fork your analysis creating a parallel output directory. 

Instead of repeating the whole analysis, you can "fork" it at the coverage step with this script.

The output directory of this script can be submitted to the original pipeline, which will begin 
from the window analysis step, using the new coverage files produced.



  -s	--source-directory		Directory where original analysis was run.
					The analysis has to have reached at least the coverage step.

  -c	--new-coverage			The new directory will contain subsets of the original
					DEPTH files, limited at this particular minimum coverage.
					It has to be at least the same as the original coverage.

	--just-link			Instead of creating a new set of coverage files, the program
					will just link the existing ones and all the *done files, 
					without asking for anything else. This is useful to 
					re-run the program at a different window size

  -p	--threads			Number of parallel processes.

  -d	--destination-directory		Directory where to put output files.
					A "fork" of the original directory will be created, containing
					only the necessary files to start over from the window analysis
					on. No mapping file will be copied.

'''.format(VERSION))

p = ap.ArgumentParser()
p.add_argument("-s", "--source-directory", required=True)
p.add_argument("-c", "--new-coverage", required=True, type=int)
p.add_argument("--just-link", action="store_true", default=False)
p.add_argument("-p", "--threads", type=int, default=1)
p.add_argument("-d", "--destination-directory", required=True)
args = p.parse_args()


### functions ### 

def list_depth_files(dirPath):

	allFiles = os.listdir(dirPath)
	covFiles = [ filename for filename in allFiles if filename[-5:] == "depth" ]

	return covFiles


def touch_donefiles(dest_directory):

	try:
		if os.path.exists(dest_directory) == False:
			os.mkdir(dest_directory)

		doneFiles = [
			"1_reference/build_hisat2_indexes.done",
			"1_reference/get_genome_file.done",
			"2_reads/isize_estimation/isize_estimation.done",
			"2_reads/link_read_files.done",
			"3_mapping/read_mapping.done",
			"3_mapping/read_filtering.done",
			"4_coverage/extract_coverage.done",
			]
	
		destDirs = [
			"1_reference", 
			"2_reads", 
			"2_reads/isize_estimation", 
			"3_mapping", 
			"4_coverage",
			]
	
		for dirPath in destDirs:
			destPath = dest_directory + "/" + str(dirPath)
			if os.path.exists(destPath) == False:
				os.mkdir(destPath)
	
		for donefile in doneFiles:
			open(dest_directory + "/" + donefile, "w").close()

		return True

	except:
		return False


def filter_coverage_file(source_cov_file, dest_cov_file, min_cov, outq, lock):

	# open old coverage file 
	# ...
	# read it 
	# ...
	# keep only positions at minimum coverage 

	INPUT = open(source_cov_file, "r")
	OUTPUT = open(dest_cov_file, "w")

	written = 0
	skipped = 0
	for line in INPUT:
		lst = line.rstrip("\n\r\b").split("\t")

		if int(lst[2]) >= min_cov:
			OUTPUT.write(line)
			written += 1

		else:
			skipped += 1

	INPUT.close()
	OUTPUT.close()

	lock.acquire()
	outq.put( (source_cov_file, dest_cov_file, written, skipped) )
	lock.release()


def link_file(source, destination):

	cmd = [ "ln",
		"-s",
		source, 
		destination ]
	
	if os.path.exists(destination) == False:
		subprocess.call(cmd)

	return True


def main():

	# get all the *.done files 
	# ...
	# get all the destination directories 
	# ...
	# create new directories if not existing 
	# ...
	# "touch" all *.done files in correct directories 
	# ...
	# get all coverage files in source directory 
	# ...
	# if new coverage: create filtered coverage files at new coverage
	# ...
	# if just link: link previous coverage files 
	# ...
	# write sum-up of the results 

	#######
	# begin
	#######

	sys.stderr.write("\n### BEGIN ###\n\n")


	#####################################
	# copy donefiles and create workspace
	##################################### 

	sys.stderr.write("[{0}] Generating environment\n".format(at()))
	status = touch_donefiles(args.destination_directory)

	if status == False:
		sys.stderr.write("ERROR: donefiles could not be generated.\n\n")
		sys.exit(1)

	sys.stderr.write("[{0}] ... DONE\n".format(at()))


	##################
	# copy genome file
	################## 

	sys.stderr.write("[{0}] Listing chromosome lengths\n".format(at()))

	copyfile(	"{0}/1_reference/genome.fa.lengths".format(args.source_directory), 
			"{0}/1_reference/genome.fa.lengths".format(args.destination_directory)	)

	sys.stderr.write("[{0}] ... DONE\n".format(at()))


	####################
	# coverage filtering
	#################### 

	if args.just_link == True:

		sys.stderr.write("[{0}] Linking source coverage files\n".format(at()))

		source_cov_dir = "{0}/4_coverage".format(args.source_directory)
		covFiles = list_depth_files(source_cov_dir)

		for covFile in covFiles:
			lst = covFile.split(".")
			(subgenome, region, coverage) = (lst[0], lst[1], int(lst[2].rstrip("x")))
			cwd = os.getcwd()
			source_cov_file = cwd + "/" + "/".join([args.source_directory, "4_coverage", covFile])
			dest_cov_filename = ".".join([subgenome, region, str(args.new_coverage) + "x", "depth"])
			dest_cov_file = cwd + "/" + "/".join([args.destination_directory, "4_coverage", dest_cov_filename])

			status = link_file(source_cov_file, dest_cov_file)

			if status == False:
				break

		if status == False:
			sys.stderr.write("ERROR: couldn't link the coverage files\n\n")
			sys.exit(1)

		else:
			open("{0}/4_coverage/filter_coverage_file.done".format(args.destination_directory), "w").close()

	else:
		sys.stderr.write("[{0}] Filtering coverage files\n".format(at()))
	
		mp.set_start_method('spawn')
		manager = mp.Manager()
		lock = manager.Lock()
		outq = manager.Queue()
	
		pool = mp.Pool(processes=args.threads)
	
		source_cov_dir = "{0}/4_coverage".format(args.source_directory)
		covFiles = list_depth_files(source_cov_dir)
	
		for covFile in covFiles:
			lst = covFile.split(".")
			(subgenome, region, coverage) = (lst[0], lst[1], int(lst[2].rstrip("x")))
			source_cov_file = "/".join([args.source_directory, "4_coverage", covFile])
			dest_cov_filename = ".".join([subgenome, region, str(args.new_coverage) + "x", "depth"])
			dest_cov_file = "/".join([args.destination_directory, "4_coverage", dest_cov_filename])
	
			if args.new_coverage < coverage:
				sys.stderr.write("ERROR: new coverage ({0}x) seems to be lower than the one in the source directory ({1}x). Can't filter a coverage file requiring a lower coverage than it was already used.\n\n".format(args.new_coverage, coverage))
				sys.exit(1)
	
			else:
				func_args = (	source_cov_file, 
						dest_cov_file, 
						args.new_coverage,
						outq, 
						lock	)
	
				pool.apply_async(	filter_coverage_file, 
							func_args	)
	
		pool.close()
		pool.join()
	
		open("{0}/4_coverage/filter_coverage_file.done".format(args.destination_directory), "w").close()

		Results = []
		while outq.qsize() > 0:
			Results.append(outq.get())
	
		Results = sorted(Results, key=itemgetter(0))
		for lst in Results:
			sys.stderr.write("\t".join([str(x) for x in lst]) + "\n")

	sys.stderr.write("[{0}] ... DONE\n".format(at()))

	#####
	# end
	#####

	sys.stderr.write("\n### END ###\n\n")
	sys.stderr.write("--- TIME: %s seconds ---\n\n" % (time.time() - start_time))


if __name__ == "__main__":

	main()
