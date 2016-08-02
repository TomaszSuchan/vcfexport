#! /usr/bin/env python3

import vcf
import re
import argparse
import string

def parse_args():
    parser = argparse.ArgumentParser(
            description='Converts VCF files to other formats, outputs to STDOUT')
    parser.add_argument('-i',
            type=str,
            dest='input',
            help='Input file.')
    parser.add_argument('-o',
            type=argparse.FileType('w'),
            dest="output",
            help='Output file.')
    parser.add_argument('-t',
            type=str,
            dest='type',
            choices=set(('fRs','nexus')),
            help='Output file type; supports nexus, and fineRADstructure (fRs) so far.')
    return parser.parse_args()

def nexus(input, output):
    IUPAC = {'AC': 'M',
             'CA': 'M',
             'AG': 'R',
             'GA': 'M',
             'AT': 'W',
             'TA': 'M',
             'CG': 'S',
             'GC': 'M',
             'CT': 'Y',
             'TC': 'M',
             'GT': 'K',
             'TG': 'M',
            }
    #list of samples
    vcf_reader = vcf.Reader(open(input, 'r'))
    sample_list = vcf_reader.samples
    sample_n = len(sample_list)
    #list for loci for each sample, each locus is [allele,allele]
    sample_concat = ['' for x in range(sample_n)]
    #this produces concatenation of all snps per sample:
    for line in vcf_reader:
        sample_nr=0
        for sample in line.samples:
            genotype = sample.gt_bases
            if genotype != None:
                base1 = re.split('\W+', genotype)[0]
                base2 = re.split('\W+', genotype)[1]
                if base1==base2:
                    sample_concat[sample_nr] += base1
                else:
                    sample_concat[sample_nr] += IUPAC[base1+base2]
            else:
                sample_concat[sample_nr] += 'N'
            sample_nr += 1
    #print sample name and contatenated seq with IUPAC codes for heterozygotes
    output.write('#NEXUS\n')
    output.write('BEGIN DATA;\n')
    output.write('  DIMENSIONS NTAX=%d NCHAR=%d;\n' % (sample_n, len(sample_concat[0])) )
    output.write('  FORMAT DATATYPE=DNA MISSING=N GAP=- INTERLEAVE=NO;\n')
    output.write('  MATRIX\n')
    max_width = len(max(sample_list, key=len))
    sample_nr = 0
    for sample in sample_concat:
        output.write('  ' + sample_list[sample_nr].ljust(max_width+2) + sample + '\n')
        sample_nr += 1
    output.write('    ;\n')
    output.write('END;')
      
def fineRADstructure(input, output):
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
                contigs[contig_nr][sample_nr][0] += ''
                contigs[contig_nr][sample_nr][1] += '' 
            sample_nr += 1
    #write output:
    output.write('\t'.join(map(str,sample_list)) + '\n') #samplenames
    for key in contigs.keys():
        contig_list = []
        for sample in contigs[key]:
            if sample[0] == '' and sample[1] == '':
                contig_list.append('')
            else:
                contig_list.append('/'.join(sample))
        output.write('\t'.join(map(str,contig_list)) + '\n')

def main():
    args = parse_args()
    if args.type == 'nexus':
        nexus(args.input, args.output)
    elif args.type == 'fRs':
        fineRADstructure(args.input, args.output)

if __name__ == "__main__":
    main()
