# MalariaART

Utility function to translate nucleotide variants to mutated protein sequences from [MalariaGen](https://www.malariagen.net/resource/34)


Example usage:

'''
# target_gene is the standard gene name such as PF3D7_1343700
db = GetVariantForGene(target_sample_ids=sample_ids)
output_dict = db.get_variant_in_gene(target_gene)
genotype_df = db.get_protein_sequence_df(output_dict)
'''

The genotype_df file will be a pandas dataframe. Each row is a mutated protein sequence for the sample for the target_gene
