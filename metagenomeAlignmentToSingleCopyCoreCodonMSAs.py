# Rauf Salamzade
# Kalan Lab, MMI, UW-Madison


"""
This code is largely adapted from our methods for lsaBGC-DiscoVary and explained in the methods for Salamzade, Swaney, and Kalan 2022. 
"""

import os
import sys
from time import sleep
import argparse
from Bio import SeqIO
from collections import defaultdict
from operator import itemgetter
import subprocess
import multiprocessing 
import pysam
import statistics
from scipy import stats
import itertools
import logging 

valid_alleles = set(['A', 'C', 'G', 'T'])

def runAlignmentApproach():
    """
    Void function which runs primary workflow for program.
    """

    """
    PARSE REQUIRED INPUTS
    """

    paired_end_sequencing_file = 'Byrd2017_Metagenomes_Test.txt' # A A file listing the location of the forward and reverse Byrd et al. 2017 metagenomic readsets.
    codon_alignments_file =  'Codon_Alignments_367.txt'   # A file listing the homolog group identifier and path to the codon alignment for it.
    outdir = os.path.abspath('DiscoVary_Results_367_Test/') + '/' # The resulting directory
    cpus = 20
    
    ### vet input files quickly
    try:
        assert (os.path.isfile(paired_end_sequencing_file))
        assert (os.path.isfile(codon_alignments_file))
    except:
        raise RuntimeError('One or more of the input files provided, does not exist. Exiting now ...')

    if os.path.isdir(outdir):
        sys.stderr.write("Output directory exists. Overwriting in 5 seconds ...\n ")
        #sleep(5)
    else:
        os.system('mkdir %s' % outdir)
    
    genes_representative_fasta = outdir + 'GCF_Genes_Representatives.fasta'
    bowtie2_db_prefix = outdir + 'GCF_Genes_Representatives'
    ref_al_to_seq = extractGenesAndCluster(genes_representative_fasta, codon_alignments_file, bowtie2_db_prefix)
   
    bowtie2_outdir = outdir + 'Bowtie2_Alignments/'
    if not os.path.isfile(bowtie2_outdir): os.system('mkdir %s' % bowtie2_outdir)
    logObject = createLoggerObject(outdir + 'Bowtie2_Log.txt')
    runBowtie2Alignments(bowtie2_db_prefix, paired_end_sequencing_file, bowtie2_outdir, logObject, cpus=cpus)

    results_outdir = outdir + 'Alignment_Parsing/'
    if not os.path.isdir(results_outdir): os.system('mkdir %s' % results_outdir)
    runSNVMining(paired_end_sequencing_file, ref_al_to_seq, genes_representative_fasta, codon_alignments_file, bowtie2_outdir, results_outdir, cpus=cpus)

    sys.exit(0)

def determineAllelesFromCodonAlignment(codon_alignment, max_mismatch=10, matching_percentage_cutoff=0.99, filter_by_genelength=True):
	gene_sequences = {}
	gene_sequences_lengths = []
	allele_identifiers = {}
	seqs_comprehensive = set([])
	with open(codon_alignment) as oca:
		for i, rec in enumerate(SeqIO.parse(oca, 'fasta')):
			gene_sequences_lengths.append(len(str(rec.seq).upper().replace('N', '').replace('-', '')))
	median_length = statistics.median(gene_sequences_lengths)
	mad_length = max(stats.median_abs_deviation(gene_sequences_lengths, scale="normal"), 5)
	with open(codon_alignment) as oca:
		for i, rec in enumerate(SeqIO.parse(oca, 'fasta')):
			gene_seq_len = len(str(rec.seq).upper().replace('N', '').replace('-', ''))
			if filter_by_genelength and (gene_seq_len < (median_length-mad_length) or gene_seq_len > (median_length+mad_length)):
				continue
			gene_sequences[rec.id] = str(rec.seq).upper()
			allele_identifiers[rec.id] = i
			seqs_comprehensive.add(rec.id)

	pairs = []
	seqs_paired = set([])
	pair_matching = defaultdict(lambda: defaultdict(float))
	for i, g1 in enumerate(gene_sequences):
		g1s = gene_sequences[g1]
		for j, g2 in enumerate(gene_sequences):
			if i >= j: continue
			g2s = gene_sequences[g2]
			tot_comp_pos = 0
			g1_comp_pos = 0
			g2_comp_pos = 0
			match_pos = 0
			mismatch_pos = 0
			for pos, g1a in enumerate(g1s):
				g2a = g2s[pos]
				if g1a in valid_alleles and g2a in valid_alleles:
					if g1a != g2a:
						mismatch_pos += 1
				if g1a in valid_alleles or g2a in valid_alleles:
					tot_comp_pos += 1
					if g1a == g2a:
						match_pos += 1
				if g1a in valid_alleles:
					g1_comp_pos += 1
				if g2a in valid_alleles:
					g2_comp_pos += 1
			general_matching_percentage = float(match_pos)/float(tot_comp_pos)
			g1_matching_percentage = float(match_pos)/float(g1_comp_pos)
			g2_matching_percentage = float(match_pos)/float(g2_comp_pos)
			pair_matching[g1][g2] = general_matching_percentage
			if general_matching_percentage >= matching_percentage_cutoff or g1_matching_percentage >= matching_percentage_cutoff or g2_matching_percentage >= matching_percentage_cutoff:
				if mismatch_pos <= max_mismatch:
					seqs_paired.add(g1)
					seqs_paired.add(g2)
					pairs.append(sorted([g1, g2]))

	"""	
	Solution for single-linkage clustering taken from mimomu's repsonse in the stackoverflow page:
	https://stackoverflow.com/questions/4842613/merge-lists-that-share-common-elements?lq=1
	"""
	L = pairs
	LL = set(itertools.chain.from_iterable(L))
	for each in LL:
		components = [x for x in L if each in x]
		for i in components:
			L.remove(i)
		L += [list(set(itertools.chain.from_iterable(components)))]

	for seq in seqs_comprehensive:
		if not seq in seqs_paired:
			L.append([seq])

	allele_cluster_min_id = {}
	for allele_cluster in L:
		gene_identifiers = set([])
		for gene in allele_cluster:
			gene_identifiers.add(allele_identifiers[gene])
		min_gi = min(gene_identifiers)
		allele_cluster_min_id[min_gi] = allele_cluster

	allele_clusters = defaultdict(set)
	for i, aci in enumerate(sorted(allele_cluster_min_id.keys())):
		for gene in allele_cluster_min_id[aci]:
			allele_clusters['Allele_Cluster_' + str(i+1)].add(gene)

	return [allele_clusters, pair_matching]


def extractGenesAndCluster(genes_representative_fasta, codon_alignments_file, bowtie2_db_prefix):
    grf_handle = open(genes_representative_fasta, 'w')
    ref_al_to_seq = defaultdict(dict)
    with open(codon_alignments_file) as ocaf:
        for line in ocaf:
            line = line.strip()
            hg, cod_alignment_fasta = line.split('\t')
            gene_to_seq = {}
            with open(cod_alignment_fasta) as ocaf:
                for rec in SeqIO.parse(ocaf, 'fasta'):
                    gene_to_seq[rec.id] = str(rec.seq).replace('-', '').replace('N', '')
            alleles_clustered, pair_matching = determineAllelesFromCodonAlignment(cod_alignment_fasta)
            for allele_cluster in alleles_clustered:
                best_rep_score = defaultdict(float)
                for ag1 in alleles_clustered[allele_cluster]:
                    for ag2 in alleles_clustered[allele_cluster]:
                        best_rep_score[ag1] += pair_matching[ag1][ag2]
                representative_gene = sorted([x for x in best_rep_score if best_rep_score[x] == max(best_rep_score.values())])[0]
                ref_al_to_seq[hg][hg + '|' + allele_cluster + '|' + representative_gene] = str(gene_to_seq[representative_gene])
                grf_handle.write('>' + hg + '|' + allele_cluster + '|' + representative_gene + '\n' + str(gene_to_seq[representative_gene]) + '\n')
    grf_handle.close()

    bowtie2_build = ['bowtie2-build', genes_representative_fasta, bowtie2_db_prefix]
    subprocess.call(' '.join(bowtie2_build), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
								executable='/bin/bash')
    return (ref_al_to_seq)

def runSNVMining(paired_end_sequencing_file, ref_al_to_seq, bowtie2_ref_fasta, codon_alignments_file, bowtie2_outdir, results_dir, cpus=1):
    gene_pos_to_msa_pos = defaultdict(lambda: defaultdict(dict))
    gene_pos_to_allele = defaultdict(lambda: defaultdict(dict))
    codon_alignment_lengths = defaultdict(int)
    with open(codon_alignments_file) as ocaf:
        for line in ocaf:
            line = line.strip()
            hg, cod_alignment = line.split('\t')
            seq_count = 0
            with open(cod_alignment) as oca:
                for j, rec in enumerate(SeqIO.parse(oca, 'fasta')):
                    gene_id = rec.id
                    real_pos = 1
                    seq_without_gaps = len([1 for bp in str(rec.seq) if bp != '-']) 
                    for msa_pos, bp in enumerate(str(rec.seq)):
                        if j == 0: codon_alignment_lengths[hg] += 1
                        if bp != '-':
                            gene_pos_to_allele[hg][gene_id][real_pos] = bp.upper()
                            gene_pos_to_msa_pos[hg][gene_id][real_pos] = msa_pos+1
                            real_pos += 1
                    seq_count += 1

    process_args = []
    with open(paired_end_sequencing_file) as opesf:
        for line in opesf:
            line = line.strip()
            sample = line.split('\t')[0]
            process_args.append([sample, bowtie2_outdir + sample + '.sorted.bam',
                                bowtie2_ref_fasta, results_dir, dict(ref_al_to_seq), 
                                dict(gene_pos_to_msa_pos), dict(gene_pos_to_allele),
                                dict(codon_alignment_lengths)])

    p = multiprocessing.Pool(cpus)
    p.map(snv_miner_single, process_args)
    p.close()

def snv_miner_single(input_args):
    sample, bam_alignment, ref_fasta, res_dir, ref_al_to_seq, gene_pos_to_msa_pos, gene_pos_to_allele, codon_alignment_lengths = input_args
    
    if not os.path.isfile(bam_alignment): return
    result_file = res_dir + sample + '.txt'
    res_outf = open(result_file, 'w')
    res_outf.write('Contig,Position,Sample-A,Sample-C,Sample-G,Sample-T\n')

    bam_handle = pysam.AlignmentFile(bam_alignment, 'rb')
    for hg in ref_al_to_seq:
        hg_ref_genes = ref_al_to_seq[hg].keys()
        read_ascpus_per_allele = defaultdict(list)
        hg_genes_covered = 0
        gene_sequence = {}
        total_reads = set([])
        for rg in hg_ref_genes: 
            gene_seq = ref_al_to_seq[hg][rg]
            gene_sequence[rg] = gene_seq
            gstart = 1
            gend = len(gene_seq)
            
            gene_length = gend - gstart + 1
            gene_covered_1 = 0

            try:
                for pileupcolumn in bam_handle.pileup(contig=rg, stepper="nofilter"):
                    pos_depth = 0
                    for pileupread in pileupcolumn.pileups:
                        if pileupread.is_del or pileupread.is_refskip: continue
                        read = pileupread.alignment
                        if read.query_qualities[pileupread.query_position] < 30: continue
                        pos_depth += 1
                        if pos_depth >= 1:
                            gene_covered_1 += 1
            except:
                pass

            gene_coverage_1 = gene_covered_1 / float(gene_length)
            if gene_coverage_1 < 0.9: continue
            hg_genes_covered += 1
            
            for read_alignment in bam_handle.fetch(rg): 
                read_name = read_alignment.query_name
                total_reads.add(read_name)
                read_ascore = read_alignment.tags[0][1]
                read_ref_positions = set(read_alignment.get_reference_positions())
                first_real_alignment_pos = None
                last_real_alignment_pos = None
                indel_positions = set([])
                matches = set([])
                for b in read_alignment.get_aligned_pairs(with_seq=True):
                    if not (b[0] == None or b[1] == None):
                        if first_real_alignment_pos == None:
                            first_real_alignment_pos = b[1]
                        last_real_alignment_pos = b[1]
                        if b[2].isupper():
                            matches.add(b[1])
                    else:
                        indel_positions.add(b[1])

                main_alignment_positions = set(range(first_real_alignment_pos, last_real_alignment_pos + 1))
                sum_indel_len = len(main_alignment_positions.intersection(indel_positions))
                matching_percentage = float(len(matches))/float(len(main_alignment_positions))

                read_ascpus_per_allele[read_name].append([rg, read_ascore, matching_percentage, len(main_alignment_positions), sum_indel_len, read_alignment])

        accounted_reads = set([])
        hg_align_pos_alleles = defaultdict(lambda: defaultdict(set))
        for read in read_ascpus_per_allele:
            top_score = -1000000
            score_sorted_alignments = sorted(read_ascpus_per_allele[read], key=itemgetter(1), reverse=True)
            for i, align in enumerate(score_sorted_alignments):
                if i == 0: top_score = align[1]
                if align[1] == top_score and ((align[2] >= 0.99 and align[3] >= 60) or (align[2] >= 0.95 and align[3] >= 100)) and align[4] <= 5:
                        read_alignment = align[-1]
                        min_read_ref_pos = min(read_alignment.get_reference_positions())
                        read_referseq = read_alignment.get_reference_sequence().upper()
                        read_queryseq = read_alignment.query_sequence
                        read_queryqua = read_alignment.query_qualities
                            
                        for b in read_alignment.get_aligned_pairs(with_seq=True):
                            if b[0] == None or b[1] == None: continue
                            ref_pos = b[1]+1
                            alt_al = read_queryseq[b[0]].upper()
                            ref_al = read_referseq[b[1] - min_read_ref_pos].upper()
                            assert (ref_al == str(gene_sequence[align[0]]).upper()[b[1]])
                            if b[2] == 'n' or ref_al == 'N' or alt_al == 'N': continue
                            que_qual = read_queryqua[b[0]]
                            if (que_qual >= 30) and ((ref_pos+3) < len(gene_sequence[align[0]])):
                                cod_pos = gene_pos_to_msa_pos[hg][align[0].split('|')[-1]][ref_pos]
                                hg_align_pos_alleles[cod_pos][alt_al].add(read)
                                accounted_reads.add(read)

        for pos in range(1, codon_alignment_lengths[hg]+1):
            printlist = [hg, str(pos)]
            for al in ['A', 'C', 'G', 'T']:
                printlist.append(str(len(hg_align_pos_alleles[pos][al])))
            res_outf.write(','.join(printlist) + '\n')

    res_outf.close()


def runBowtie2Alignments(bowtie2_reference, paired_end_sequencing_file, bowtie2_outdir, logObject, cpus=1):
	"""
	Wrapper function for running Bowtie2 alignments to reference database/index
	:param bowtie2_reference: path to Bowtie2 reference/index
	:param paired_end_sequencing_file: tab delimited file with three columns: (1) sample name (2) path to forward
									   reads and (3) path to reverse reads
	:param bowtie2_outdir: Path to directory where Bowtie2 results should be written
	:param logObject: logging object for documentation
	:param cpus: number of cpus (total) to use. If more than 4 cpus provided, then parallel Bowtie2 jobs with 4 cpus
				  each will be started.
	"""
	bowtie2_cpus = cpus
	bowtie2_pool_size = 1
	if cpus >= 4:
		bowtie2_cpus = 4
		bowtie2_pool_size = int(cpus / 4)

	try:
		bowtie2_inputs = []
		with open(paired_end_sequencing_file) as opesf:
			for line in opesf:
				line = line.strip()
				sample = line.split('\t')[0]
				reads = line.split('\t')[1:]
				bowtie2_inputs.append([sample, reads, bowtie2_reference, bowtie2_outdir, bowtie2_cpus, logObject])
		p = multiprocessing.Pool(bowtie2_pool_size)
		p.map(bowtie2_alignment, bowtie2_inputs)
		p.close()
	except Exception as e:
		logObject.error("Issues in setting up and running Bowtie2 alignments.")
		logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())

def bowtie2_alignment(input_args):
	"""
	Function to perform Bowtie2 alignment of paired-end reads to a database/reference and post-processing of alignment
	file with samtools (e.g. convert to BAM format, sort BAM file, and index it).
	"""
	sample, reads, bowtie2_reference, bowtie2_outdir, bowtie2_cpus, logObject = input_args

	sam_file = bowtie2_outdir + sample + '.sam'
	bam_file = bowtie2_outdir + sample + '.bam'
	bam_file_sorted = bowtie2_outdir + sample + '.sorted.bam'
	#bam_file_filtered = bowtie2_outdir + sample + '.filtered.bam'
	#am_file_filtered_sorted = bowtie2_outdir + sample + '.filtered.sorted.bam'

	bowtie2_cmd = ['bowtie2', '--very-sensitive-local', '--no-unal', '-a', '-x', bowtie2_reference, '-U',
				   ','.join(reads), '-S', sam_file, '-p', str(bowtie2_cpus)]

	samtools_view_cmd = ['samtools', 'view', '-h', '-Sb', sam_file, '>', bam_file]
	samtools_sort_cmd = ['samtools', 'sort', '-@', str(bowtie2_cpus), bam_file, '-o', bam_file_sorted]
	samtools_index_cmd = ['samtools', 'index', bam_file_sorted]

	try:
		run_cmd(bowtie2_cmd, logObject)
		run_cmd(samtools_view_cmd, logObject)
		run_cmd(samtools_sort_cmd, logObject)
		run_cmd(samtools_index_cmd, logObject)

		os.system("rm -f %s %s" % (sam_file, bam_file)) # bam_file_sorted, bam_file_filtered, bam_file_sorted + '.bai'))
	except Exception as e:
		if bowtie2_outdir != "" and sample != "":
			os.system('rm -f %s/%s*' % (bowtie2_outdir, sample))
		raise RuntimeError(traceback.format_exc())

def run_cmd(cmd, logObject, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL):
	"""
	Simple function to run a single command through subprocess with logging.
	"""
	logObject.info('Running the following command: %s' % ' '.join(cmd))
	try:
		subprocess.call(' '.join(cmd), shell=True, stdout=stdout, stderr=stderr, executable='/bin/bash')
		logObject.info('Successfully ran: %s' % ' '.join(cmd))
	except Exception as e:
		logObject.error('Had an issue running: %s' % ' '.join(cmd))
		logObject.error(traceback.format_exc())
		raise RuntimeError('Had an issue running: %s' % ' '.join(cmd))

def createLoggerObject(log_file):
	"""
	Function which creates logger object.
	:param log_file: path to log file.
	:return: logging logger object.
	"""
	logger = logging.getLogger('task_logger')
	logger.setLevel(logging.DEBUG)
	# create file handler which logs even debug messages
	fh = logging.FileHandler(log_file)
	fh.setLevel(logging.DEBUG)
	# create console handler with a higher log level
	ch = logging.StreamHandler()
	ch.setLevel(logging.DEBUG)
	# create formatter and add it to the handlers
	formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', "%Y-%m-%d %H:%M")
	fh.setFormatter(formatter)
	ch.setFormatter(formatter)
	# add the handlers to the logger
	logger.addHandler(fh)
	logger.addHandler(ch)
	logger.addHandler(logging.StreamHandler())

	return logger

if __name__ == '__main__':
    runAlignmentApproach()
