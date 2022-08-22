import os
import sys
from collections import defaultdict
from Bio import SeqIO

proteomes_dir = '../Proteomes/'
cod_listing_file = 'Codon_Alignments_1250.txt'

kefir_prots = set([])
tuber_prots = set([])
aurim_prots = set([])
for f in os.listdir(proteomes_dir):
    if 'kefirresidentii' in f:
        with open(proteomes_dir + f) as opf:
            for rec in SeqIO.parse(opf, 'fasta'):
                kefir_prots.add(rec.id)
    elif 'aurim' in f:
        with open(proteomes_dir + f) as opf:
            for rec in SeqIO.parse(opf, 'fasta'):
                aurim_prots.add(rec.id)
    elif 'tubercul' in f:
        with open(proteomes_dir + f) as opf:
            for rec in SeqIO.parse(opf, 'fasta'):
                tuber_prots.add(rec.id)

top367 = set([x.strip() for x in open('Top_367_Most_Discriminatory_SCC_HGs.txt').readlines()])

with open(cod_listing_file) as oclf:
    for line in oclf:
        line = line.strip()
        hg, cod_align = line.split('\t')
        if not hg in top367: continue
        
        kefir_seqs = []
        tuber_seqs = []
        aurim_seqs = []
        ckefir_seqs = []
        ctuber_seqs = []
        caurim_seqs = []
        with open(cod_align) as oca:
            for rec in SeqIO.parse(oca, 'fasta'):
                if rec.id in kefir_prots:
                    kefir_seqs.append(list(str(rec.seq)))
                else:
                    ckefir_seqs.append(list(str(rec.seq)))
                if rec.id in tuber_prots:
                    tuber_seqs.append(list(str(rec.seq)))
                else:
                    ctuber_seqs.append(list(str(rec.seq)))
                if rec.id in aurim_prots:
                    aurim_seqs.append(list(str(rec.seq)))
                else:
                    caurim_seqs.append(list(str(rec.seq)))

        kefir_scc_sites = {}
        for i, als in enumerate(zip(*kefir_seqs)):
            als_set = set(list(als))
            if len(als_set) == 1 and not '-' in als_set:
                kefir_scc_sites[i] = list(als)[0]

        non_unique = set([])
        for i, als in enumerate(zip(*ckefir_seqs)):
            for al in set(list(als)):
                if i in kefir_scc_sites and al  == kefir_scc_sites[i]:
                    non_unique.add(i)

        for pos in kefir_scc_sites:
            if not pos in non_unique:
                print('kefir\t' + hg + '\t' + str(pos+1) + '\t' + kefir_scc_sites[pos])

        tuber_scc_sites = {}
        for i, als in enumerate(zip(*tuber_seqs)):
            als_set = set(list(als))
            if len(als_set) == 1 and not '-' in als_set:
                tuber_scc_sites[i] = list(als)[0]

        non_unique = set([])
        for i, als in enumerate(zip(*ctuber_seqs)):
            for al in set(list(als)):
                if i in tuber_scc_sites and al == tuber_scc_sites[i]:
                    non_unique.add(i)

        for pos in tuber_scc_sites:
            if not pos in non_unique:
                print('tuber\t' + hg + '\t' + str(pos+1) + '\t' + tuber_scc_sites[pos])

        aurim_scc_sites = {}
        for i, als in enumerate(zip(*aurim_seqs)):
            als_set = set(list(als))
            if len(als_set) == 1 and not '-' in als_set:
                aurim_scc_sites[i] = list(als)[0]

        non_unique = set([])
        for i, als in enumerate(zip(*caurim_seqs)):
            for al in set(list(als)):
                if i in aurim_scc_sites and al == aurim_scc_sites[i]:
                    non_unique.add(i)

        for pos in aurim_scc_sites:
            if not pos in non_unique:
                print('aurim\t' + hg + '\t' + str(pos+1) + '\t' + aurim_scc_sites[pos])




