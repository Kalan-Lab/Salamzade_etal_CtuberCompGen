import os
import sys
from collections import defaultdict
import statistics

align_dir = 'DiscoVary_Results_367/Alignment_Parsing_Looser/'
meta_file = 'Metagenomic_Profiling/Byrd_et_al_2017_MetaData.txt'
species_sites = 'Species_Specific_Alleles.txt'
scalar_file = 'Metagenome_Base_Normalization_Scalars.txt'

scalars = {}
for line in open(scalar_file):
    line = line.strip()
    ls = line.split('\t')
    scalars[ls[0]] = ls[1]

grouping = {}
scorad = {}
with open(meta_file) as omf:
    for line in omf:
        line = line.strip()
        ls = line.split('\t')
        group = ls[4]
        samp = ls[0]
        grouping[samp] = group
        scorad[samp] = '0.0'
        if ls[5] != 'NA':
            scorad[samp] = ls[5]

total_species_sites = defaultdict(int)
kefir_alleles = defaultdict(dict)
tuber_alleles = defaultdict(dict)
aurim_alleles = defaultdict(dict)
with open(species_sites) as oksf:
    for line in oksf:
        line = line.strip()
        sp, hg, pos, al = line.split('\t')
        total_species_sites[sp] += 1
        if sp == 'kefir':
            kefir_alleles[hg][int(pos)] = al.upper()
        elif sp == 'tuber':
            tuber_alleles[hg][int(pos)] = al.upper()
        elif sp == 'aurim':
            aurim_alleles[hg][int(pos)] = al.upper()

print('\t'.join(['sample', 'group', 'scorad', 'scalar', 'kefir_mc', 'tuber_mc', 'aurim_mc', 'statistic']))
for f in os.listdir(align_dir):
    m = f.split('.txt')[0]
    kefir_site_coverages = []
    tuber_site_coverages = []
    aurim_site_coverages = []
    with open(align_dir + f) as oaf:
        for i, line in enumerate(oaf):
            if i == 0: continue
            line = line.strip()
            hg, pos, a, c, g, t = line.split(',')
            pos = int(pos)
            if pos in kefir_alleles[hg]:
                if kefir_alleles[hg][pos] == 'A':
                    kefir_site_coverages.append(int(a))
                elif kefir_alleles[hg][pos] == 'C':
                    kefir_site_coverages.append(int(c))
                elif kefir_alleles[hg][pos] == 'G':
                    kefir_site_coverages.append(int(g))
                elif kefir_alleles[hg][pos] == 'T':
                    kefir_site_coverages.append(int(t))
            if pos in tuber_alleles[hg]:
                if tuber_alleles[hg][pos] == 'A':
                    tuber_site_coverages.append(int(a))
                elif tuber_alleles[hg][pos] == 'C':
                    tuber_site_coverages.append(int(c))
                elif tuber_alleles[hg][pos] == 'G':
                    tuber_site_coverages.append(int(g))
                elif tuber_alleles[hg][pos] == 'T':
                    tuber_site_coverages.append(int(t))
            if pos in aurim_alleles[hg]:
                if aurim_alleles[hg][pos] == 'A':
                    aurim_site_coverages.append(int(a))
                elif aurim_alleles[hg][pos] == 'C':
                    aurim_site_coverages.append(int(c))
                elif aurim_alleles[hg][pos] == 'G':
                    aurim_site_coverages.append(int(g))
                elif aurim_alleles[hg][pos] == 'T':
                    aurim_site_coverages.append(int(t))
   
    kefir_mc = statistics.mean(kefir_site_coverages)
    tuber_mc = statistics.mean(tuber_site_coverages)
    aurim_mc = statistics.mean(aurim_site_coverages)
    ratio_stat = (kefir_mc + 1)/(tuber_mc + aurim_mc + 1)
    print(m + '\t' + grouping[m] + '\t' + scorad[m] + '\t' + scalars[m] + '\t' + str(kefir_mc) + '\t' + str(tuber_mc) + '\t' + str(aurim_mc) + '\t' + str(ratio_stat))
