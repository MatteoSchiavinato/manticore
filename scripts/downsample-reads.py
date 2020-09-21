#!/usr/bin/env python3

import argparse as ap
from Bio import SeqIO
import subprocess 
from math import ceil 
from time import asctime as at 
import sys 

# help
if len(sys.argv) == 1:
	sys.argv.append("--help")

if sys.argv[1] in ["-h", "--help", "-help", "getopt", "usage"]:
	sys.exit('''

  --reference		FASTA file to use to compute coverage				[-]
  --read-1		READ_1 file (first in pair)					[-]
  --read-2		READ_2 file (second in pair)					[-]
  --target-coverage	Reads will be downsampled to this coverage 			[mandatory]
  --outdir		Downsampled reads will be put in this directory 		[current]
  --basename		Downsampled read files will have this prefix 			[ds]
  --seqtk-path		Path to the seqtk executable (if not in `$PATH`)		[seqtk]
  --seqtk-seed		Seed for random reads selection (same seed: same selection)	[15]


''')

# parser
p = ap.ArgumentParser()
p.add_argument("--reference", help="FASTA file to compute coverage from", required=True)
p.add_argument("--read-1", help="Read 1 FASTQ file", required=True)
p.add_argument("--read-2", help="Read2 FASTQ file")
p.add_argument("--target-coverage", help="Coverage to downsample the reads to [INT, example: 150]", type=int, required=True)
p.add_argument("--outdir", help="Output directory", default=".")
p.add_argument("--basename", help="Prefix for output files", default="ds")
p.add_argument("--seqtk-path", default="seqtk")
p.add_argument("--seqtk-seed", default=15)
args = p.parse_args()


### functions ### 

def run_cmd(cmd, stdout_file, stderr_file):

	STDOUT_FILE = open(stdout_file, "w")
	STDERR_FILE = open(stderr_file, "w")
	code = subprocess.call(cmd, stdout=STDOUT_FILE, stderr=STDERR_FILE)
	STDOUT_FILE.close()
	STDERR_FILE.close()

	return code


def calculate_genome_size(reference):

	# compute total length 
	INPUT = open(reference, "r")
	Lengths = [ len(str(record.seq)) for record in SeqIO.parse(INPUT, "fasta") ]
	total_length = sum(Lengths)
	del Lengths
	INPUT.close()

	return int(total_length)


def compute_read_total_length(reads_file):

	# compute total length
	INPUT = open(reads_file, "r")
	Lengths = [ len(str(record.seq)) for record in SeqIO.parse(INPUT, "fastq") ]
	total_length = sum(Lengths)
	total_num = len(Lengths)
	del Lengths
	INPUT.close()

	return (int(total_length), int(total_num))


def downsample_reads(read_file_1, read_file_2, ratio, outdir, basename, target_cov, seqtk_path, seed):

	Read_files = {
		read_file_1 : "{0}/{1}.{2}x.1.fastq".format(outdir, basename, target_cov),
		read_file_2 : "{0}/{1}.{2}x.2.fastq".format(outdir, basename, target_cov)
	}

	for read_file in Read_files.keys():

		cmd = [	"seqtk", 
			"sample", 
			"-s", str(seed),
			str(read_file),
			str(ratio) ]

		code = run_cmd(	cmd, 
					"{0}".format(Read_files[read_file]), 
					"{0}.stderr".format(Read_files[read_file]) )

		if code != 0:
			break

	if code == 0:
		return True
	else:
		return False


# main script 
if __name__ == "__main__":

	# compute genome size
	genome_size = calculate_genome_size(args.reference)
	sys.stderr.write("[{0}] Genome size computed: {1} Gbp\n".format(at(), genome_size/1e9))

	# compute read total length 
	# (save total number to use later)
	(read_1_length, read_1_num) = compute_read_total_length(args.read_1)
	(read_2_length, read_2_num) = compute_read_total_length(args.read_2)
	sys.stderr.write("[{0}] READ 1: {1} reads, total length {2} nt\n".format(at(), 	read_1_num, 
											read_1_length ))
	sys.stderr.write("[{0}] READ 2: {1} reads, total length {2} nt\n".format(at(), 	read_2_num, 
											read_2_length ))

	if read_1_num != read_2_num:
		sys.exit("ERROR: you have a different number of read 1 and 2; something is wrong.\n\n")

	# compute current coverage 
	total_length = read_1_length + read_2_length
	current_cov = float(total_length) / float(genome_size)
	sys.stderr.write("[{0}] Reads currently cover the provided genome with ~{1}x coverage\n".format(at(), round(current_cov, 2)))

	# compute ratio desired / current 
	ratio = float(args.target_coverage) / float(current_cov)
	sys.stderr.write("[{0}] {1} of the total reads will be randomly chosen to downsample the files\n".format(at(), round(ratio, 2)))
	if ratio > float(1):
		sys.exit("ERROR: you are asking for more coverage ({0}) than these reads can afford ({1}x).\n\n".format(args.target_coverage, round(current_cov, 2)))

	# re-parse reads and write to output those with idx 
	status = downsample_reads(	args.read_1, 
					args.read_2, 
					ratio, 
					args.outdir, 
					args.basename, 
					args.target_coverage,
					args.seqtk_path,
					args.seqtk_seed )
	if status == True:
		sys.stderr.write("[{0}] Reads were downsampled to {1}x coverage.\n\n".format(at(), args.target_coverage))
	elif status == False:
		sys.stderr.write("ERROR: reads could not be downsampled\n\n")
		sys.exit(1)
