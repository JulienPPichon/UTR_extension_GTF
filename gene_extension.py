"""This script extend the 3'UTRs region of each transcript to have a better mapping of library with a high number of reads in 3'UTR regions."""

import argparse
import copy
import re
import sys
from datetime import datetime
from tqdm import tqdm
from intervaltree import Interval, IntervalTree

def get_args(argv = None):
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--filename", help="take filename.gtf as input")
	parser.add_argument("-o", "--output", help="name of the gtf output file with extended genes", default = default_name_output)
	parser.add_argument("-l", "--extent_length", type=int, help="Length of last exon transcript extension", default = 2000)
	parser.add_argument("-d", "--gene_distance", type=int, help="minimum distance separing 2 transcripts of different genes", default = 50)
	parser.add_argument("-e", "--exon_number", help="add exon number at the end of the attributes", action="store_true")
	parser.add_argument("-c", "--chromosome_size", help="Use chromosome size file to extend last gene", default = "no_input_file")
	return parser.parse_args(argv)
	

def gtf_to_dict(filename):
	"""Add exon lines to a dictionnary with chromosome number as key and intervaltree object.
	
	Convert each line into a list with "tab" and ";" as delimiters. If the third column mentions an exon, add it to the dict_exons dictionnary,
	otherwise to the dict_other one.
	If the first column is not in the key of the dictionnary, add it and convert it to an intervaltree. If there is, add the value as interval 
	object (attributes as a list of each column) to the corresponding key, leading to a dictionnary of list of interval.
	Check and append to a list exons with a value of the start higher than the end.
	"""
	dict_exons = {}
	# Dictionnary containning chromosome as key and intervaltree of list of exons as values. 
	dict_all = {}
	# Junk dictionnary containning lines other that exon type.  
	header = []
	# Save the header in case there is one in the gtf.
	exons_wrong_interval = []  
	# List of exons for which its end is before its start.
	with open(filename, "r") as filin:
		for line in filin:
			if line.startswith("#!"):
				header.append(line)
			else:			
				line = line.strip()
				line = re.split("\t|;", line)[:-1]
				# Transform the string line in a list type with tab and ; as delimiters. 
				if line[0] not in dict_exons.keys():
				# Check if the chromosome is already in the dictionnary.  
					dict_exons[str(line[0])] = []
					dict_exons[str(line[0])] = IntervalTree()
					dict_all[str(line[0])] = []
				if int(line[3]) < int(line[4]): 
				# Check if the end position of the exon is after the start one.
					if line[2] == "exon":
						dict_exons[str(line[0])][int(line[3]):int(line[4])] = line
					dict_all[str(line[0])].append(line)
				else:
					exons_wrong_interval.append(line)
		dict_exons = copy.deepcopy(dict_exons)
	return [dict_all, dict_exons, header, exons_wrong_interval]


def dict_chromosome_size(chromosome_size):
	"""Create a dictionnary of chromosomes sizes thanks to -c input."""
	dict_size = {}
	with open(chromosome_size, "r") as filin_size:
		for line in filin_size:
			line = line.strip()
			line = re.split("\t", line)
			dict_size[str(line[0])] = line[1]
	return dict_size
	

def find_pos(dict_exons):
	"""Determine the position of gene_id and transcript_id in the list of attributes of an exon."""
	for n, attributes in enumerate(sorted(dict_exons[list(dict_exons.keys())[0]])[0].data):
		if attributes.startswith(" gene_id ") or attributes.startswith("gene_id "):
			gene_col = n
		if attributes.startswith(" transcript_id ") or attributes.startswith("transcript_id "):
			transcript_col = n
	return [gene_col, transcript_col]


def list_to_set_wi(exons_wrong_interval):
	genes_wrong_interval = set()
	for exons in exons_wrong_interval:
		genes_wrong_interval.add(exons[gene_col])
	return genes_wrong_interval	
	

def calculate_interval(dict_exons, gene_col, transcript_col):
	"""If the requirements are respected, extend the last exons end of each transcripts.
	
	Iterate on each exons of each chromosomes with the start position sorted.
	For each loop, keep in memory the strand and gene_id of the previous exons.
	If the strand of the current gene is positive, add the exon number, making a double loop from the start of
	the gene to look at the transcript id.
	If the gene_id change, look back to the previous exon to summon exon number for negative strand,
	and determine if extension of the lasts exons for each transcript is possible for the previous gene.
	Contrary to the positive strand, negative strand exon assignation of a gene is performed after determined
	the end of the gene, because the last exon is the one with the smaller position (i.e. already iterate).
	A reverse loop is done to deal with that issue.
	Different requirement meet 4 possibilies for extension. If previous and current genes are + strand, add the input
	length extension if the last exon of the previous gene are still spaced of gene_distance input value from the first exon 
	from the current gene; otherwise, add the maximum length allowed, preventing the gene overlappings and keeping the space.
	If the previous gene is + and the current one is -, extension risks to overlapped. So, the requirements check the double
	of the gene extension that should be performed. If requirements are not respected, extend by the half length between the 2 genes
	minus the half allowed gene_distance input.
	The same kind of operation are made in the "- +" and "- -" cases.
	Two kind of errors due to gtf format are often meet in the step.
	Add to a list the exon for which the start begin before the end of the previous gene.
	Add to a list the exon if its gene_id has already been meet previously and be interrupted by an other gene.
	"""
	overlapping_genes = set()
	exon_trans = 1
	duplicated_genes = set()	
	already_known = set()
	warning_genes = set()
	nb_genes = set()
	# When the gene_id change, add to this set the gene_id of the previous exon. 
	for chromosome in tqdm(dict_exons):	
		data_exon = sorted(dict_exons[chromosome])[0].data
		# First exon attributes of the chromosome.
		nb_genes.add(data_exon[gene_col])
		gene_strand = sorted(dict_exons[chromosome])[0].data[6]
		gene_start = sorted(dict_exons[chromosome])[0].begin
		exon_trans = 1
		gene_end = 0
		for current_exon in sorted(dict_exons[chromosome]):	
			if current_exon.data[gene_col] == data_exon[gene_col]:
			# Start the case if the following exon is from the same gene and + strand.
				if args.exon_number:
				# Start the optionnal argument of exon number summons for + strand.
					if current_exon.data[6] == "+":			
						for exon_in_gene in sorted(dict_exons[chromosome][gene_start:current_exon.begin]):
							if current_exon.data[transcript_col] == exon_in_gene.data[transcript_col]:
								exon_trans += 1
						current_exon.data.append("exon : " + str(exon_trans))
						exon_trans = 1
						# End the exon number summons.
			if current_exon.data[gene_col] != data_exon[gene_col]:
				nb_genes.add(current_exon.data[gene_col])
			# Start the case if the following exon is not from the same gene than the previous one.
				if current_exon.data[gene_col] in already_known:
					duplicated_genes.add(data_exon[gene_col])
					warning_genes.add(data_exon[gene_col])
					continue
					# Add the exon to a list if it has already been see previously.
				if gene_start >= current_exon.begin or sorted(dict_exons[chromosome][gene_start:current_exon.begin])[-1].end >= current_exon.begin:
					overlapping_genes.add(current_exon.data[gene_col])
					overlapping_genes.add(data_exon[gene_col])
					warning_genes.add(current_exon.data[gene_col])
					warning_genes.add(data_exon[gene_col])			
					continue
					# Skip the exon if it's start begin before the end of the previous gene.
				prv_gene_end = copy.deepcopy(gene_end)
				# Keep in memory the gene_end of the previous gene.
				prv_gene_strand = copy.deepcopy(gene_strand)
				# Keep in memory the strand of the previous gene.
				gene_strand = sorted(dict_exons[chromosome][gene_start:current_exon.begin])[-1].data[6]	
				gene_end = sorted(dict_exons[chromosome][gene_start:current_exon.begin])[-1].end
				# Determine the end and the strand of the gene to extend.
				if args.exon_number:
				# Start the optionnal argument of exon number summons for - strand.
					if data_exon[6] == "-":
						exon_trans = 0
						for exon_in_gene in sorted(dict_exons[chromosome][gene_start:gene_end + 1], reverse = True):
							for following_exons in sorted(dict_exons[chromosome][exon_in_gene.begin:gene_end + 1], reverse = True):
								if exon_in_gene.data[transcript_col] == following_exons.data[transcript_col]:
									exon_trans += 1
							exon_in_gene.data.append("exon : " + str(exon_trans))
							exon_trans = 0
					if current_exon.data[6] == "+":
					# Case where the gene change but it's + strand.
						exon_trans = 1
						current_exon.data.append("exon : " + str(exon_trans))						
				# Start the gene extension.
				if data_exon[6] == "+" and data_exon[gene_col] not in warning_genes:
					if current_exon.data[6] == "+":
					# Case where both previous and current gene strand are +.
						for exon_in_gene in sorted(dict_exons[chromosome][gene_start:gene_end + 1]):
							for past_exon in sorted(dict_exons[chromosome][gene_start:exon_in_gene.begin]):
								if past_exon.data[transcript_col] == exon_in_gene.data[transcript_col]:
									break
									# Go to the next exon of the gene if this transcript are already been extended. 
							else:
								for following_exons in sorted(dict_exons[chromosome][exon_in_gene.begin:gene_end + 1]):
									if exon_in_gene.data[transcript_col] == following_exons.data[transcript_col]:
										prv_gene_end = following_exons
										# Determine the last exon of the transcript.
								distance_gene = current_exon.begin - prv_gene_end.end
								if distance_gene > args.extent_length + args.gene_distance:
								# If not overlapping because extension.
									for reverse_sorted_exon in sorted(dict_exons[chromosome][gene_start:gene_end + 1], reverse = True):
										if reverse_sorted_exon.data[transcript_col] == exon_in_gene.data[transcript_col]:
											reverse_sorted_exon.data[4] = str(int(reverse_sorted_exon.end) + args.extent_length)
											break
											# Extend the last exon of the corresponding transcript.
								else:
								# If overlapping because extension, calculate the max extension allowed.
									for reverse_sorted_exon in sorted(dict_exons[chromosome][gene_start:gene_end + 1], reverse = True):
										if reverse_sorted_exon.data[transcript_col] == exon_in_gene.data[transcript_col] and distance_gene > args.gene_distance:
											reverse_sorted_exon.data[4] = str(int(reverse_sorted_exon.end) + distance_gene - args.gene_distance)
											break				
					if current_exon.data[6] == "-":
					# Case where current gene are - strand and previous one + strand.
						for exon_in_gene in sorted(dict_exons[chromosome][gene_start:gene_end + 1]):
							for past_exon in sorted(dict_exons[chromosome][gene_start:exon_in_gene.begin]):
								if past_exon.data[transcript_col] == exon_in_gene.data[transcript_col]:
									break
							else:
								for following_exons in sorted(dict_exons[chromosome][exon_in_gene.begin:gene_end + 1]):
									if exon_in_gene.data[transcript_col] == following_exons.data[transcript_col]:
										prv_gene_end = following_exons
								distance_gene = current_exon.begin - prv_gene_end.end
								if distance_gene > (args.extent_length * 2) + args.gene_distance:			 
									for reverse_sorted_exon in sorted(dict_exons[chromosome][gene_start:gene_end + 1], reverse = True):
										if reverse_sorted_exon.data[transcript_col] == exon_in_gene.data[transcript_col]:
											reverse_sorted_exon.data[4] = str(int(reverse_sorted_exon.end) + args.extent_length)
											break
								else:
									for reverse_sorted_exon in sorted(dict_exons[chromosome][gene_start:gene_end + 1], reverse = True):
										if reverse_sorted_exon.data[transcript_col] == exon_in_gene.data[transcript_col] and distance_gene > args.gene_distance:
											reverse_sorted_exon.data[4] = str(int(reverse_sorted_exon.end) + int(distance_gene / 2) - int(args.gene_distance / 2))
											break			
				if data_exon[6] == "-" and data_exon[gene_col] not in warning_genes:
				# Case where previous gene are - strand.
					for exon_in_gene in sorted(dict_exons[chromosome][gene_start:gene_end + 1]):
						for past_exon in sorted(dict_exons[chromosome][gene_start:exon_in_gene.begin]):
							if past_exon.data[transcript_col] == exon_in_gene.data[transcript_col]:
								break
						else:
							if exon_in_gene.begin == dict_exons[chromosome].begin():
								if exon_in_gene.begin > args.extent_length:
									exon_in_gene.data[3] = str(int(exon_in_gene.begin) - args.extent_length)
								else:
									exon_in_gene.data[3] = str(1)
							else:
								if prv_gene_strand == "+":
								# Case where current gene are + strand.
									if int(prv_gene_end) < exon_in_gene.begin - (args.extent_length * 2) - int(args.gene_distance / 2):
										exon_in_gene.data[3] = str(int(exon_in_gene.begin) - args.extent_length)
									else:
										if int(exon_in_gene.begin - int(prv_gene_end)) > args.gene_distance:
											exon_in_gene.data[3] = str((int(exon_in_gene.begin) - int((exon_in_gene.begin - int(prv_gene_end)) / 2)) + int(args.gene_distance / 2))
								if prv_gene_strand == "-":
								# Case where current gere are - strand.
									if int(prv_gene_end) < exon_in_gene.begin - args.extent_length - args.gene_distance:
					 					exon_in_gene.data[3] = str(int(exon_in_gene.begin) - args.extent_length)
									else:
										if int(exon_in_gene.begin - int(prv_gene_end)) > args.gene_distance:
											exon_in_gene.data[3] = str((int(exon_in_gene.begin) - int((exon_in_gene.begin - int(prv_gene_end)) )) + args.gene_distance )		
				# End the gene extension.			
				already_known.add(data_exon[gene_col])
				data_exon = current_exon.data[:-1]
				gene_start = current_exon.begin
				#dict_transcript[current_exon.data[transcript_col]] = [gene_start]			
			if current_exon.end == dict_exons[chromosome].end():
			# Case where the exon is the last one of the chromosome.
				nb_genes.add(current_exon.data[gene_col])
				prv_gene_end = copy.deepcopy(gene_end)
				last_gene_end = current_exon.end
				if current_exon.data[6] == "+":
					for exon__in_gene in sorted(dict_exons[chromosome][gene_start:last_gene_end + 1]):
						for past_exon in sorted(dict_exons[chromosome][gene_start:exon_in_gene.begin]):
							if past_exon.data[gene_col] == exon_in_gene.data[gene_col]:
								continue
						for reverse_sorted_exon in sorted(dict_exons[chromosome][gene_start:last_gene_end + 1], reverse = True):
							if args.chromosome_size == "no_input_file":
								if reverse_sorted_exon.data[transcript_col] == exon_in_gene.data[transcript_col]:
									reverse_sorted_exon.data[4] = str(int(reverse_sorted_exon.end) + args.extent_length)
									break
							else:
								if reverse_sorted_exon.data[transcript_col] == exon_in_gene.data[transcript_col]:
									if int(reverse_sorted_exon.end) + args.extent_length < int(dict_size[current_exon.data[0]]):
										reverse_sorted_exon.data[4] = str(int(reverse_sorted_exon.end) + args.extent_length)
									else:
										reverse_sorted_exon.data[4] = str(dict_size[current_exon.data[0]])
									break																
				if current_exon.data[6] == "-":
					if args.exon_number:
						exon_trans = 0
						for reverse_sorted_exon in sorted(dict_exons[chromosome][gene_start:last_gene_end + 1], reverse = True):
							for reverse_sorted_exon2 in sorted(dict_exons[chromosome][reverse_sorted_exon.begin:last_gene_end + 1], reverse = True):
								if reverse_sorted_exon.data[transcript_col] == reverse_sorted_exon2.data[transcript_col]:
									exon_trans += 1
							reverse_sorted_exon.data.append("exon : " + str(exon_trans))
							exon_trans = 0
					for exon_in_gene in sorted(dict_exons[chromosome][gene_start:last_gene_end + 1]):
						for past_exon in sorted(dict_exons[chromosome][gene_start:exon_in_gene.begin]):
							if past_exon.data[transcript_col] == exon_in_gene.data[transcript_col]:
								break
						else:
							if exon_in_gene.begin == dict_exons[chromosome].begin():
								continue
							else:
								if prv_gene_strand == "+":
									if int(prv_gene_end) < exon_in_gene.begin - (args.extent_length * 2) - (args.gene_distance / 2):
										exon_in_gene.data[3] = str(int(exon_in_gene.begin) - args.extent_length)
									else:
										if int(exon_in_gene.begin) - int(prv_gene_end) > args.gene_distance:
											exon_in_gene.data[3] = str(int(int(exon_in_gene.begin) - int((exon_in_gene.begin - int(prv_gene_end)) / 2)) + int(args.gene_distance / 2) )
								if prv_gene_strand == "-":
									if int(prv_gene_end) < exon_in_gene.begin - args.extent_length - args.gene_distance:
				 						exon_in_gene.data[3] = str(int(exon_in_gene.begin) - args.extent_length)
									else:
										if int(exon_in_gene.begin) - int(prv_gene_end) > args.gene_distance:
											exon_in_gene.data[3] = str(int(int(exon_in_gene.begin) - int((exon_in_gene.begin - int(prv_gene_end)) )) + args.gene_distance )
				already_known.add(current_exon.data[gene_col])
	return [overlapping_genes, duplicated_genes, dict_exons, warning_genes, already_known, nb_genes]


def merge_dict(dict_all, dict_exons, gene_col, transcript_col, warning_genes):
	"""Change the exons of the original gtf with the extended ones without warnings."""
	cpt = 0
	for chromosome in dict_all:
		gene_start = 1
		gene_end = 1000000
		for line in dict_all[chromosome]:
			if line[2] == "exon" and line[gene_col] not in warning_genes:
				previous_line = line
				break
		for line in dict_all[chromosome]:
			if line[2] == "exon" and line[gene_col] not in warning_genes:
				for exon in sorted(dict_exons[chromosome][int(line[3]):int(line[4]) + 1]):
					if line[6] == "+":
						if args.exon_number:
							if exon.data[:4] == line[:4] and exon.data[5:-1] == line[5:]:
								line[4] = exon.data[4]
								line.append(exon.data[-1])
						else:
							if exon.data[:4] == line[:4] and exon.data[5:] == line[5:]:
								line[4] = exon.data[4]
					if line[6] == "-":
						if args.exon_number:
							if exon.data[:3] == line[:3] and exon.data[4:-1] == line[4:]:
								line[3] = exon.data[3]
								line.append(exon.data[-1])
						else:
							if exon.data[:3] == line[:3] and exon.data[4:] == line[4:]:
								line[3] = exon.data[3]
				previous_line = copy.deepcopy(line)
	return dict_all
	

def dict_to_gtf(dict_all):
	"""Create a output file similar to the gtf input with the gene extended."""
	with open(args.output, "w") as filout:
		for line in header:
			filout.write("{}".format(line))
		for chromosome in dict_all:		
			for exon in dict_all[chromosome]:
				exon_gtf = "\t".join(exon[0:9]) + "; " + "; ".join(exon[9:]) + ";"
				filout.write("{}\n".format(exon_gtf))


def errors_to_resume(genes_wrong_interval, overlapping_genes, duplicated_genes):
	"""Create an output file resuming the encountered errors during the script processing."""
	name_resume = sys.argv[2].rsplit('.',1)[0] + ".extent.txt"
	if args.output:
		name_resume = args.output + ".txt"
	with open(name_resume, "w") as filout2:
		filout2.write("{} genes with at least an exon with start finishing after the end (these exons have been removed from the gtf):\n\n".format(len(genes_wrong_interval)))
		for gene in genes_wrong_interval:
			filout2.write("{}, ".format(gene))
		filout2.write("\n\n{} genes overlapping another (i.e. end of a gene finishing after the start of the next gene, extension not performed.):\n\n".format(len(overlapping_genes)))
		for gene in overlapping_genes:
			filout2.write("{}, ".format(gene))
		filout2.write("\n\n{} genes for which the gene_id has already been see and been interrupted by an other gene_id \
		(extension not performed). Most of the time, this is also a case of overlapping genes:\n\n".format(len(duplicated_genes)))
		for gene in duplicated_genes:
			filout2.write("{}, ".format(gene))


if __name__ == "__main__":

	
	print(str(datetime.now()))
	if len(sys.argv) == 1:
		default_name_output = "filename.extent.gtf"
	if len(sys.argv) >= 2:
		if sys.argv[1] == "-h":
			default_name_output = "filename.extent.gtf"
		else:
			default_name_output = sys.argv[2].rsplit('.',1)[0] + ".extent.gtf"
	else:
		default_name_output = "filename.extent.gtf"
	argvals = None
	args = get_args(argvals)	
	dict_all, dict_exons, header, exons_wrong_interval = gtf_to_dict(args.filename)
	if args.chromosome_size != "no_input_file":
		dict_size = dict_chromosome_size(args.chromosome_size)
	gene_col, transcript_col = find_pos(dict_exons)
	genes_wrong_interval = list_to_set_wi(exons_wrong_interval)
	overlapping_genes, duplicated_genes, dict_exons, warning_genes, already_known, nb_genes = calculate_interval(dict_exons, gene_col, transcript_col)
	for i in genes_wrong_interval:
		warning_genes.add(i)
	for y in genes_wrong_interval:
		nb_genes.add(y)
	dict_all = merge_dict(dict_all, dict_exons, gene_col, transcript_col, warning_genes)
	dict_to_gtf(dict_all)		
	errors_to_resume(genes_wrong_interval, overlapping_genes, duplicated_genes)	
	warning_number = len(genes_wrong_interval) + len(overlapping_genes) + len(duplicated_genes)
	print(str(datetime.now()))
	if warning_number != 0:
		print("{} warnings, please check txt file.\n{} of {} genes have been totally or partially extended.".format(warning_number, len(nb_genes) - len(warning_genes), len(nb_genes)))

	
