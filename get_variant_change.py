from collections import OrderedDict
import pickle

import re
import numpy as np
import pandas as pd

from Bio import SeqIO
from Bio.Data import IUPACData
from Bio.Seq import Seq


import malariagen_data
import os

# from proteinbert import OutputType, OutputSpec, FinetuningModelGenerator, load_pretrained_model, finetune, evaluate_by_len
# from proteinbert.conv_and_global_attention_model import get_model_with_hidden_layers_as_outputs

class GetVariantForGene():

	one_to_three = IUPACData.protein_letters_1to3
	three_to_one = IUPACData.protein_letters_3to1

	def __init__(self, 
			     target_sample_ids):
		"""
		target_sample_ids: a list of sample ids, this is necessary to reduce the size of the xarray_genotype. At its minimum, all samples didn't pass QC should be removed.
		"""
		# get CDS file for Pf3D7
		self.protein_db = self._get_protein_db()
		# pf7 api from malariagen
		self.pf7 = malariagen_data.Pf7()
		self.variant_meta_df, self.xarray_genotype = self._parse_pf7()
		self.target_sample_ids = target_sample_ids

		self.sample_ids = self.xarray_genotype.sample_id.values
		self.keep_sample_idx = pd.Series(self.sample_ids).isin(self.target_sample_ids)
		self.keep_sample_ids = self.sample_ids[self.keep_sample_idx]
		self.variant_calls = self.pf7.variant_calls(extended=True)

	def get_variant_in_gene(self, target_gene):
		pickle_file = f'{target_gene}_{len(self.target_sample_ids)}.pkl'
		# if the variant file exsits, load it
		if os.path.isfile(pickle_file):
			with open(pickle_file, 'rb') as output_file:
				output_dict = pickle.load(output_file)
				return output_dict
		# otherwise, generate the variant file
		else:
			gene_info = self.protein_db[target_gene]
			target_chrom = gene_info['chrom']
			target_start = gene_info['start']
			target_end = gene_info['end']
			strand = gene_info['strand']
			ref_seq = gene_info['seq']
			idx = ((self.variant_meta_df['variant_chrom'] == target_chrom) & 
				(self.variant_meta_df['variant_position'] >= target_start) & 
				(self.variant_meta_df['variant_position'] <= target_end))

			target_genotype = self.xarray_genotype[np.where(idx)[0][::6]//6, :] # becasue idx is a boolean array corresponding to df_amino, which has 6 rows for each variant so we need to divide by 6
			target_aa = self.variant_meta_df[idx]
			# output_dict = {i:{'sample_id':j, 'NA_change':[], 'AA_change':[]} for i,j in enumerate(self.target_sample_ids)}
			output_dict = {i:{'sample_id':j, 'NA_change':[]} for i,j in enumerate(self.target_sample_ids)} # NA_change stores nucleotide change
			output_dict.update({'ref_seq': ref_seq})

			for i, (variants, _) in enumerate(target_aa.iterrows()):
				alt_idx = i % 6
				if alt_idx != 0:
					continue

				target_aa_chunk = target_aa.iloc[i:i+6]

				genotype_data = self.variant_calls.call_genotype[variants[0],self.keep_sample_idx].values
				unique_genotypes = np.unique(genotype_data)

				# skip if all genotypes are -1 or 0
				if (unique_genotypes < 1).all():
					continue

				# target_aa_change = target_aa_chunk['variant_ANN_HGVS_p'].map(self._get_amino_acid_change).dropna().tolist()
				target_na_change = target_aa_chunk['variant_ANN_HGVS_c'].map(self._get_nucleic_acid_change).dropna().tolist()
				if strand == '+':
					for counter, v in enumerate(target_na_change):
						sample_idx = np.where(np.max(genotype_data, axis=1)==(counter+1))[0]
						for sample_idx_single in sample_idx:
							# output_dict[sample_idx_single]['AA_change'].append(v)
							output_dict[sample_idx_single]['NA_change'].append(v)
				else:
					for counter, v in enumerate(reversed(target_na_change)):
						sample_idx = np.where(np.max(genotype_data, axis=1)==(counter+1))[0]
						for sample_idx_single in sample_idx:
							# output_dict[sample_idx_single]['AA_change'].append(v)
							output_dict[sample_idx_single]['NA_change'].append(v)

			with open(f'{target_gene}_{len(self.target_sample_ids)}.pkl', 'wb') as output_file:
				pickle.dump(output_dict, output_file)

			return output_dict
		
	def _get_amino_acid_change(self, variant_ANN_HGVS_p):
		match = re.findall('([A-Z][a-z]+)([0-9]+)([A-Z][a-z]+)', variant_ANN_HGVS_p)
		if match == []:
			return None
		else:
			aa_ref, aa_pos, aa_alt = match[0] 
			return self.three_to_one.get(aa_ref, None), aa_pos, self.three_to_one.get(aa_alt, None)

	def _get_nucleic_acid_change(self, variant_ANN_HGVS_c):
		match = re.findall('([0-9]+)([A-Z]+)>([A-Z]+)', variant_ANN_HGVS_c)
		if match == []:
			return None
		else:
			nt_pos, nt_ref, nt_alt = match[0]
			return nt_pos, nt_ref, nt_alt

	def _get_nucleic_acid_change_reverse(self, variant_ANN_HGVS_c):
		MAP = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
		match = re.findall('([0-9]+)([A-Z]+)>([A-Z]+)', variant_ANN_HGVS_c)
		if match == []:
			return None
		else:
			nt_pos, nt_ref, nt_alt = match[0] 
			return nt_pos, MAP[nt_ref], MAP[nt_alt]
				
	def get_protein_sequence_df(self, output_dict):
		sequence_df = {}
		ref_seq = output_dict['ref_seq']
		for k in range(len(self.target_sample_ids)):
			v = output_dict[k]
			alt_seq = list(ref_seq)
			for na_change in v['NA_change']:
				na_pos, na_ref, na_alt = na_change
				assert ref_seq[int(na_pos)-1] == na_ref, "Reference nucleotide does not match!"
				alt_seq[int(na_pos)-1] = na_alt
			seq = self._translate_cds_to_protein(''.join(alt_seq))
			sequence_df[v['sample_id']] = seq if seq[-1]!='*' else seq[:-1]

		return pd.DataFrame.from_dict(sequence_df, orient='index')

	def _parse_pf7(self):
		meta = self.pf7.sample_metadata()
		keep_samples = meta['Sample'][meta['QC pass']==True].values
		population = meta['Population'][meta['QC pass']==True]
		vc = self.pf7.variant_calls(extended=True)
		alt_alleles = vc.variant_ANN_HGVS_c.to_dataframe()
		alt_alleles_p = vc.variant_ANN_HGVS_p.to_dataframe()
		ann_df = vc.variant_ANN_Annotation.to_dataframe()
		vc_cat = ann_df['variant_ANN_Annotation'].astype('category')
		variants_dict = OrderedDict()
		variants_dict['Type'] = pd.Series(vc['variant_is_snp'][:]).map({True: 'SNP', False: 'non-SNP'})
		variants_dict['Multiallelic'] = pd.Series(vc['variant_numalt'][:] == 1).map({True: 'Bi-allelic', False: 'Multi-allelic'})
		variants_dict['Coding'] = pd.Series(vc['variant_CDS'][:]).map({True: 'Coding', False: 'Non-coding'})
		variants_dict['is_pass'] = vc['variant_filter_pass'][:]
		variants_dict['num_alleles'] = vc['variant_numalt'][:] + 1
		variants_dict['Chrom'] = pd.Series(vc['variant_chrom'][:])
		variants_dict['Pos'] = pd.Series(vc['variant_position'][:])
		df_variants = pd.DataFrame.from_dict(variants_dict)
		df_variants.loc[df_variants['Type']=='non-SNP', 'Multiallelic'] = ''
		keep_idx = (df_variants['is_pass']==True) & (df_variants['Type']=='SNP') & (df_variants['Coding']=='Coding')
		df = df_variants[keep_idx]
		df_amino = alt_alleles_p[alt_alleles_p.index.get_level_values('variants').isin(df.index)]
		df_nucleotide = alt_alleles[alt_alleles.index.get_level_values('variants').isin(df.index)]
		df_amino['variant_ANN_HGVS_c'] = df_nucleotide['variant_ANN_HGVS_c'].copy()

		xarray_genotype = vc.call_genotype[df.index]
		return df_amino, xarray_genotype

	def _translate_cds_to_protein(self, cds):
		coding_dna = Seq(cds)
		protein_seq = coding_dna.translate()
		return str(protein_seq)

	# def _get_alt_alleles_p(chrom, start, end):
	# 	return alt_alleles_p[(alt_alleles_p['variant_chrom'] == chrom) & (alt_alleles_p['variant_position'] >= start) & (alt_alleles_p['variant_position'] <= end)]


	def _get_protein_db(self, use_protein_annot=False):
		# release 66
		if use_protein_annot:
			fasta_file = "/media/changli/Share3/MalariaGen/malariaART/data/PlasmoDB-66_Pfalciparum3D7_AnnotatedProteins.fasta"  # protein
		else:
			fasta_file = "/media/changli/Share3/MalariaGen/malariaART/data/PlasmoDB-66_Pfalciparum3D7_AnnotatedCDSs.fasta"  # cds
		# Parse the fasta file
		protein_db = {}
		for record in SeqIO.parse(fasta_file, "fasta"):
			fasta_description = record.description
			# Split the description by '|' to get all fields and trim spaces
			fields = [x.strip() for x in fasta_description.split('|')]
			# Now loop over the fields to find location
			for field in fields:
				if 'location' in field:
					location_info = field.split('=')[1]
				if 'SO' in field:
					so_info = field.split('=')[1]
			gene_name = fields[0].split('.')[0]
			if so_info == 'protein_coding_gene':
				chrom, pos = location_info.split(':')
				start, end, strand, *_ = re.split(r'[-()]', pos)
				protein_db[gene_name] = {'chrom': chrom, 'start':int(start), 'end':int(end), 'seq':str(record.seq), 'strand':'-' if strand=='' else strand}
		
		return protein_db