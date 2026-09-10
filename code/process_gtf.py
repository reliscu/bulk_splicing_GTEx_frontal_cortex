import numpy as np
import pandas as pd
from gtfparse import read_gtf

def process_gtf(gtf_file, exclude, no_trim_id, gene_type_tag, transcript_type_tag):

    print('Processing GTF file...')
    
    exclude_chromosomes = exclude.split(',')
    
    gtf = read_gtf(gtf_file)
    gtf = gtf.to_pandas()

    try:

        try:
            cols = ['seqname', 'start', 'end', 'feature', 'strand', 'transcript_id', gene_type_tag, 'gene_name', 'gene_id', transcript_type_tag, 'exon_id', 'exon_number', 'frame', 'tag', 'protein_id']
            gtf = gtf.loc[pd.Series([x not in exclude_chromosomes for x in gtf.seqname]), cols]
            gtf.columns = ['chrom', 'start', 'end', 'feature', 'strand', 'transcript', 'gene_type', 'gene_name', 'gene_id', 'transcript_type', 'exon_id', 'exon_number', 'frame', 'tag', 'protein_id']
        except:
            cols = ['seqname', 'start', 'end', 'feature', 'strand', 'transcript_id', 'gene_name', 'gene_id', 'exon_id', 'exon_number', 'frame', 'tag', 'protein_id']
            gtf = gtf.loc[pd.Series([x not in exclude_chromosomes for x in gtf.seqname]), cols]
            gtf.columns = ['chrom', 'start', 'end', 'feature', 'strand', 'transcript', 'gene_name', 'gene_id', 'exon_id', 'exon_number', 'frame', 'tag', 'protein_id']
            
    except:
        raise Exception('Isufficient information to create annotation. transcript_id is needed to find cassette exons.')

    if not no_trim_id:
        gtf['gene_id'] = gtf.gene_id.str.replace(r'\.\d+$', '', regex=True)

    gtf.transcript = [x.split('.')[0] for x in gtf.transcript]
    gtf.exon_id = [x.split('.')[0] for x in gtf.exon_id]
     
    gtf = gtf.loc[gtf.gene_id != '']
        
    return gtf