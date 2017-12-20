#! /usr/bin/env python3

# in base1 = re.split('\W+', genotype)[0] use | if phased!

import vcf
import re
import argparse
import string
from random import randint
from Bio import SeqIO

IUPAC = {'AC': 'M',
         'CA': 'M',
         'AG': 'R',
         'GA': 'R',
         'AT': 'W',
         'TA': 'W',
         'CG': 'S',
         'GC': 'S',
         'CT': 'Y',
         'TC': 'Y',
         'GT': 'K',
         'TG': 'K'
        }

def parse_args():
    parser = argparse.ArgumentParser(
            description='Converts VCF files to other formats.')
    parser.add_argument('-i',
            type=str,
            dest='input',
            help='Input file.')
    parser.add_argument('-r',
            type=str,
            dest='ref_input',
            help='Reference file.')
    parser.add_argument('-o',
            type=argparse.FileType('w'),
            dest="output",
            help='Output file.')
    parser.add_argument('-t',
            type=str,
            dest='type',
            choices=set(('fRs', 'structure','nexus','fasta','phylip')),
            help='Output file type; supports fine RAD structure, structure, nexus, fasta, and phylip format.')
    parser.add_argument('--all',
            dest="all",
            action="store_true",
            default=False,
            help='Output all  SNPs per contig (default: onle unlinked SNP per contig)')
    return parser.parse_args()


#def linkedSNPs(input, reference):
#    #read fasta reference
#    handle = open(reference, "rU")
#    reference_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
#    handle.close()
#    
#    #list for contigs, sequences in each loci inside
#    samples_sequences = [{key: list(reference_dict[key].seq) for key in contigs_list} for x in range(sample_n)]
#    
#    #free memory from the full reference
#    del(reference_dict)
#    
#    #this produces sequences per sample per contig:
#    vcf_reader = vcf.Reader(open(input, 'r'))
#    for line in vcf_reader:
#        contig_nr = line.CHROM
#        sample_nr=0
#        for sample in line.samples:
#            genotype = sample.gt_bases
#            position = sample.site.POS
#            if genotype != None:
#                base1 = re.split('\W+', genotype)[0]
#                base2 = re.split('\W+', genotype)[1]
#                if base1==base2:
#                    SNP = base1
#                else:
#                    SNP = IUPAC[base1+base2]
#            else:
#                SNP = 'N'
#            #position starts from 1:
#            #TODO: what to do with more complex variants?
#            samples_sequences[sample_nr][contig_nr][position+1] = SNP
#            sample_nr += 1
#    #join sequences fom list to strings here:
#    samples_sequences = [{key: ''.join(sample[key]) for key in sample} for sample in samples_sequences]
#    
#    return([samples_sequences, sample_list])

def SNPs(input, concat, contigs_list, sample_n):
    #list for contigs, concat SNPs in each loci inside
    samples = [{key: [] for key in contigs_list} for x in range(sample_n)]
    #this produces concatenation of all snps per sample per contig, phased:
    vcf_reader = vcf.Reader(open(input, 'r'))
    for line in vcf_reader:
        contig_nr = line.CHROM
        #contig_nr = contig_nr.split('_')[2]
        #contig_nr = int(contig_nr)
        sample_nr=0
        for sample in line.samples:
            genotype = sample.gt_bases
            if genotype != None:
                base1 = re.split('\W+', genotype)[0]
                base2 = re.split('\W+', genotype)[1]
                if base1=='N' or base2=='N':
                    if base1!='N':
                        samples[sample_nr][contig_nr] += base1
                    elif base2!='N':
                        samples[sample_nr][contig_nr] += base2
                    else:
                        samples[sample_nr][contig_nr] += '-'
                elif base1==base2:
                    samples[sample_nr][contig_nr] += base1
                else:
                    samples[sample_nr][contig_nr] += IUPAC[base1+base2]
            else:
                samples[sample_nr][contig_nr] += '-'
            sample_nr += 1
    if concat==False:
        #select single, unlinked SNP
        #TODO: select SNP from full columns
        loci_lengths=[len(samples[0][d]) for d in samples[0]]
        random_snp=[randint(0,l-1) for l in loci_lengths]
        samples_concat=[]
        for sample_nr in range(sample_n):
            concatenated=''
            nr=0
            for contig in samples[sample_nr]:
                concatenated += samples[sample_nr][contig][random_snp[nr]]
                nr += 1
            samples_concat.append(concatenated)
        return(samples_concat)
    elif concat==True:
        #return concatenated SNPs:
        #loci_lengths=[len(samples[0][d]) for d in samples[0]]
        #random_snp=[randint(0,l-1) for l in loci_lengths]
        samples_concat=[]
        for sample_nr in range(sample_n):
            concatenated=''
            nr=0
            for contig in samples[sample_nr]:
                concatenated += ''.join(samples[sample_nr][contig])
                nr += 1
            samples_concat.append(concatenated)
        return(samples_concat)

def nexus(input, output, concat, contigs_list, sample_n, sample_list):
    #get concetenated unlinked or linked SNPs
    samples_concat = SNPs(input, concat, contigs_list, sample_n)
    
    #print sample name and contatenated seq in nexus format
    output.write('#NEXUS\n')
    output.write('BEGIN DATA;\n')
    output.write('  DIMENSIONS NTAX=%d NCHAR=%d;\n' % (sample_n, len(samples_concat[0])) )
    output.write('  FORMAT DATATYPE=DNA MISSING=N GAP=- INTERLEAVE=NO;\n')
    output.write('  MATRIX\n')
    max_width = len(max(sample_list, key=len))
    sample_nr = 0
    for sample in samples_concat:
        output.write('  ' + sample_list[sample_nr].ljust(max_width+2) + sample + '\n')
        sample_nr += 1
    output.write('    ;\n')
    output.write('END;')


def fasta(input, output, concat, contigs_list, sample_n, sample_list):
    samples_concat = SNPs(input, concat, contigs_list, sample_n)
    
    #print sample name and contatenated seq in fasta format
    sample_nr = 0
    for sample in samples_concat:
        output.write('>' + sample_list[sample_nr] + '\n')
        output.write(sample + '\n')
        sample_nr += 1

def structure(input, output, type):
    vcf_reader = vcf.Reader(open(input, 'r'))
    #list of samples
    sample_list = vcf_reader.samples
    sample_n = len(sample_list)
    #list of contigs
    contigs_list = []
    for line in vcf_reader:
        contig_nr = line.CHROM
        contig_nr = contig_nr.split('_')[2]
        contig_nr = int(contig_nr)
        contigs_list.append(contig_nr)
    contigs_list = list(set(contigs_list))
    contigs_n = len(contigs_list)
    # print stats
    print('Processing file with', contigs_n, 'contigs and', sample_n, 'individuals...')
    #list for contigs, concat SNPs in each loci inside
    contigs = {key: [['',''] for x in range(sample_n)] for key in contigs_list}
    #this produces concatenation of all snps per sample per contig, phased:
    vcf_reader = vcf.Reader(open(input, 'r'))
    for line in vcf_reader:
        contig_nr = line.CHROM
        contig_nr = contig_nr.split('_')[2]
        contig_nr = int(contig_nr)
        sample_nr=0
        for sample in line.samples:
            genotype = sample.gt_bases
            if genotype != None:
                contigs[contig_nr][sample_nr][0] += re.split('\W+', genotype)[0]
                contigs[contig_nr][sample_nr][1] += re.split('\W+', genotype)[1]
            else:
                contigs[contig_nr][sample_nr][0] += 'N'
                contigs[contig_nr][sample_nr][1] += 'N' 
            sample_nr += 1
    if type=='fineRADstructure':
        output.write('\t'.join(map(str,sample_list)) + '\n') #samplenames
        for key in contigs.keys():
            contig_list = []
            for sample in contigs[key]:
                contig_list.append('/'.join(sample))
            output.write('\t'.join(map(str,contig_list)) + '\n')
    elif type=='structure':
        STRUCTURE = {'A': '1',
                     'C': '2',
                     'T': '3',
                     'G': '4',
                     'N': '-9'
                    }
        # generate random SNP numbers
        # (in dictionary contig_nr:random_SNP_position):    
        random_snp = {key: [randint(0,len(contigs[key][0][0])-1)] for key in contigs.keys()}
        #print sample name and unlinked SNPs in structure format
        sample_nr = 0
        population_nr = 0
        population = ''
        for sample in sample_list:
            if population != sample_list[sample_nr].split('_')[0]:
                population_nr += 1
            population = sample_list[sample_nr].split('_')[0]
            # each sample has two lines:
            # first line:
            #    name + 5 fields according to faststructure requirements
            output.write(sample_list[sample_nr] +'\t' + str(population_nr) +'\t'+'\t'+'\t'+'\t')
            for key in contigs.keys():
                position = random_snp[key][0]
                snp = contigs[key][sample_nr][0][position]
                snp = STRUCTURE[snp]
                output.write('\t' + snp)
            output.write('\n')
            # second line:
            output.write(sample_list[sample_nr] + '\t' + str(population_nr) +'\t'+'\t'+'\t'+'\t')
            for key in contigs.keys():
                position = random_snp[key][0]
                snp = contigs[key][sample_nr][1][position]
                snp = STRUCTURE[snp]
                output.write('\t' + snp)
            output.write('\n')
            sample_nr += 1

def main():
    args = parse_args()
    #list of samples
    vcf_reader = vcf.Reader(open(args.input, 'r'))
    sample_list = vcf_reader.samples
    sample_n = len(sample_list)
    #list of contigs
    contigs_list = []
    for line in vcf_reader:
        contig_nr = line.CHROM
        #contig_nr = contig_nr.split('_')[2]
        #contig_nr = int(contig_nr)
        contigs_list.append(contig_nr)
    contigs_list = list(set(contigs_list))
    contigs_n = len(contigs_list)
    # print stats
    print('Processing file with', contigs_n, 'contigs and', sample_n, 'samples...')
    if args.type == 'nexus':
        nexus(args.input, args.output, args.all, contigs_list, sample_n, sample_list)
    elif args.type == 'fasta':
        fasta(args.input, args.output, args.all, contigs_list, sample_n, sample_list)
    elif args.type == 'structure':
        structure(args.input, args.output, type='structure')
    elif args.type == 'fRs':
        structure(args.input, args.output, type='fineRADstructure')


if __name__ == "__main__":
    main()
