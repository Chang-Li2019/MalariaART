# Author: Chang Li
# Date: 2024-02-01
# Description: Using MalariaGen data and API to Convert DNA variants and sequences to consensus protein sequences
# License: MIT License
# Copyright (c) 2024 Chang Li
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import pandas as pd
import re
from tqdm import tqdm
import os
import xarray as xr

import malariagen_data

def prepare_datafiles(protein_fa, cds_fa, genome_fa, gff_file):
    # SeqIO read fasta to dict
    protein_tmp_dict = SeqIO.to_dict(SeqIO.parse(protein_fa, "fasta"))
    protein_dict = {}
    for k, v in protein_tmp_dict.items():
        protein_dict[k.split('.')[0]] = v

    # read cds fasta
    cds_tmp_dict = SeqIO.to_dict(SeqIO.parse(cds_fa, "fasta"))
    cds_dict = {}
    for k, v in cds_tmp_dict.items():
        cds_dict[k.split('.')[0]] = v

    # read genome fasta
    genome_dict = SeqIO.to_dict(SeqIO.parse(genome_fa, "fasta"))

    # read gff
    import gffutils
    # Create a database from the GFF file
    gff_db = gffutils.create_db(gff_file, dbfn='gff.db', force=True, merge_strategy='merge', keep_order=True)
    return protein_dict, cds_dict, genome_dict, gff_db

def check_list(l):
    return all([i=='' or len(i)==1 for i in l])

def parse_desc(desc):
    d = {}
    for i in desc.split('|')[1:]:
        k,v = i.split('=')
        d[k.strip()] = v.strip()
    return d

def get_gene_seq(ref_seq, gt_col, start):
    seq = ref_seq.copy()
    for k,v in gt_col.to_dict().items():
        seq[k-start]=v
    return ''.join(seq)

def prepare_malarigen_data(species):
    if species == 'pf7':
        mg_data = malariagen_data.Pf7()
        meta = mg_data.sample_metadata()
        variant_calls = mg_data.variant_calls()
    elif species == 'ag3':
        mg_data = malariagen_data.Ag3()
        meta = mg_data.sample_metadata()
        variant_calls = mg_data.variant_calls()
    elif species == 'amin1':
        mg_data = malariagen_data.Amin1()
        meta = mg_data.sample_metadata()
        variant_calls = mg_data.variant_calls()
    elif species == 'pv4':
        mg_data = malariagen_data.Pv4()
        meta = mg_data.sample_metadata()
        variant_calls = mg_data.variant_calls()
    # quality filter
    qc_idx = variant_calls['variant_filter_pass'].values
    # snp_idx = (variant_calls['variant_is_snp']).values # too conservative by removing all * alt alleles
    snp_idx = xr.apply_ufunc(check_list, variant_calls['variant_allele'].load(), input_core_dims=[['alleles']], vectorize=True, dask='parallelized', output_dtypes=[bool]).values

    variant_allele_qc = variant_calls['variant_allele'][qc_idx&snp_idx]
    # replace * with N in variant_allele_qc
    def replace_value(entry):
        return [item.replace('*', 'N') for item in entry]
    variant_allele_qc = xr.apply_ufunc(replace_value, variant_allele_qc, input_core_dims=[['alleles']], vectorize=True, dask='parallelized', output_dtypes=[object]).values

    call_genotype_qc = variant_calls['call_genotype'][qc_idx&snp_idx]
    variant_position_qc = variant_calls['variant_position'][qc_idx&snp_idx]
    variant_chrom_qc = variant_calls['variant_chrom'][qc_idx&snp_idx]

    return meta, variant_allele_qc, call_genotype_qc, variant_position_qc, variant_chrom_qc

def argparser():
    parser = argparse.ArgumentParser(description='Convert DNA variant and sequences to protein sequences')
    # add parser argument to select from a list of species --species
    choices = ['pf7', 'ag3', 'amin1', 'pv4']
    parser.add_argument("--species", choices=choices, help="Select a species from ['pf7', 'ag3', 'amin1', 'pv4']")
    parser.add_argument('--protein_fa', type=str, default='', help='protein fasta file')
    parser.add_argument('--cds_fa', type=str, help='cds fasta file')
    parser.add_argument('--genome_fa', type=str, help='genome fasta file')
    parser.add_argument('--gff_file', type=str, help='gff file')
    parser.add_argument('--output', type=str, default='all_seqs', help='output directory')
    parser.add_argument('--log', type=str, default='log.txt', help='log file')
    args = parser.parse_args()
    return args

def main():
    # make sure all 
    args = argparser()
    protein_dict, cds_dict, genome_dict, gff_db = prepare_datafiles(args.protein_fa, args.cds_fa, args.genome_fa, args.gff_file)
    meta, variant_allele_qc, call_genotype_qc, variant_position_qc, variant_chrom_qc = prepare_malarigen_data(args.species)
    print(f'Finished database preparation for {args.species}')
    log = open(f'{args.log}', 'a')
    # if output directory exists, continue. If not, create it
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    for gene_name,info in tqdm(cds_dict.items()):
        description = parse_desc(info.description)
        chrom = re.findall(r'(.+):', description.get('location', None))[0]
        pos = re.findall(r'(\d+)-(\d+)', description.get('location', None))
        strand = re.findall(r'\([+-]\)', description.get('location', None))[0]
        length = int(description.get('length'))
        # check cds length is correct
        assert int(pos[0][1]) - int(pos[0][0]) + 1 == length
        start = int(pos[0][0])  # gene start
        end = int(pos[0][1])   # gene end
        idx = ((variant_chrom_qc==chrom) & (variant_position_qc>=start) & (variant_position_qc<=end)).values  # index for all variants of the gene, including introns
        if idx.sum() == 0:
            print(f'No variant found for {gene_name}', file=log)
            continue
        gene_va = variant_allele_qc[idx] # variant alleles
        gt_dict = {}
        for i, j in enumerate(gene_va):
            gt_dict[i] = {}
            for k, l in enumerate(j):
                if l != '':
                    gt_dict[i][k] = l 
        gene_cg = call_genotype_qc[idx].values
        # numpy array replace -1 with 0
        gene_cg[gene_cg==-1] = 0
        gene_cg_max = pd.DataFrame(gene_cg.max(axis=-1))

        allele_df = gene_cg_max.apply(lambda row: row.map(gt_dict[row.name]), axis=1)
        allele_df.index = variant_position_qc[idx].values
        ref_seq = list(genome_dict[chrom].seq[start-1:end])
        alt_seqs_dna = allele_df.apply(lambda col: get_gene_seq(ref_seq, col, start), axis=0)

        ## ref_seq may contain introns, so no need to check at this step
        # if strand == '(-)':
        #     assert Seq(''.join(ref_seq)).reverse_complement().translate()[:-1] == protein_dict[gene_name].seq
        # else:
        #     assert Seq(''.join(ref_seq)).translate()[:-1] == protein_dict[gene_name].seq

        gene_entry = gff_db[gene_name]
        mRNAs = gff_db.children(gene_entry, featuretype='mRNA')
        for mRNA in mRNAs:
            mRNA_start = mRNA.start
            mRNA_stop = mRNA.stop
            if mRNA_start == int(pos[0][0]) and mRNA_stop == int(pos[0][1]):
                canonical_mRNA = mRNA
                break

        exons = gff_db.children(canonical_mRNA, featuretype="CDS")

        exons = list(exons)
        if len(exons) == 0:
            # save the print log to a file using print function
            print(f'No CDS found for {gene_name}', file=log)
            continue
        # sort exons by start position for + and - strands, lowest first
        exon_positions = [(exon.start, exon.end) for exon in exons]
        exon_positions = sorted(exon_positions, key=lambda x: x[0])
        
        intron_coords = [(exon_positions[i][1]+1, exon_positions[i+1][0]-1) for i in range(len(exon_positions)-1)]
        # relative intron coordinates
        intron_coords = [(i[0]-exon_positions[0][0]+1, i[1]-exon_positions[0][0]+1) for i in intron_coords]
        intron_lengths = np.sum([exon_positions[i+1][0] - exon_positions[i][1] -1 for i in range(len(exon_positions)-1)])

        for i in reversed(intron_coords):
            alt_seqs_dna = alt_seqs_dna.str.slice_replace(i[0]-1, i[1], '')
            ref_seq = ref_seq[:i[0]-1] + ref_seq[i[1]:]
        if strand == '(-)':
            alt_seqs = alt_seqs_dna.map(lambda x: str(Seq(x).reverse_complement().translate())[:-1])
        else:
            alt_seqs = alt_seqs_dna.map(lambda x: str(Seq(x).translate())[:-1])
        alt_seqs.index = meta.Sample
        file_to_write = os.path.join(args.output, f'{gene_name}.csv')
        alt_seqs.to_csv(file_to_write, header=None)
        print(f'Finished {gene_name}', file=log)

    log.close()
    
if __name__ == "__main__":
    main()
