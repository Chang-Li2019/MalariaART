# MalariaART

Utility function to translate nucleotide variants to consensus protein sequences from [MalariaGen](https://www.malariagen.net/resource/34), at https://github.com/malariagen/malariagen-data-python

## Required packages:
  Biopython
  numpy
  pandas
  xarray
  malariagen_data
  
## Required arguments: 
  --species {pf7,ag3,amin1,pv4}
                        Select a species from ['pf7', 'ag3', 'amin1', 'pv4']
                        
  --protein_fa PROTEIN_FA
                        protein fasta file
                        
  --cds_fa CDS_FA       cds fasta file
  
  --genome_fa GENOME_FA
                        genome fasta file
                        
  --gff_file GFF_FILE   gff file
  
  --output OUTPUT       output directory
  
  --log LOG             log file
  
## Example usage:

```bash
python malairagen_dna_to_protein.py\
--species pf7 --log log.txt --output output_seqs --gff_file PlasmoDB-66_Pfalciparum3D7.gff \
--genome_fa PlasmoDB-66_Pfalciparum3D7_Genome.fasta --cds_fa PlasmoDB-66_Pfalciparum3D7_AnnotatedCDSs.fasta \
--protein_fa PlasmoDB-66_Pfalciparum3D7_AnnotatedProteins.fasta
```

Outputs will be saved to *output*, each file is a gene, and each row in the file is a sample.

