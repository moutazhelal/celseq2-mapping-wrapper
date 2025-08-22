#!/usr/bin/env python3

# Derived from: Mouse gastruloids transcriptomics analysis
# Upstream: https://github.com/anna-alemany/mouseGastruloids_scRNAseq_tomoseq
# Original authors (c) upstream contributors; see upstream history
# Modifications (c) 2025 Moutaz Helal
# License: GPL-3.0
# SPDX-License-Identifier: GPL-3.0-only

import sys, os
from pandas.io.parsers import read_csv
import numpy as np
import pandas as pd

try:
    gtffile = sys.argv[1]
    output = sys.argv[2]
except:
    sys.exit("Please, give: (1) gtf file; (2) prefix for output files")

#### Convert gtf file in an easy to handle data frame ####
df = read_csv(gtffile, comment='#', header=None, sep = '\t', low_memory = False)
df.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
df['attribute'] = [{f.rsplit()[0]:f.rsplit()[1].replace('"','') for f in df.loc[idx,'attribute'].rsplit(';') if len(f.rsplit())==2} for idx in df.index]

ats = set()
for idx in df.index:
    ats.update(df.loc[idx, 'attribute'].keys())

for c in ats:
    df[c] = [a[c] if c in a else '-' for a in df['attribute']]

del df['attribute']
df['gene_biotype'] = [''.join([y.capitalize() for y in x.rsplit('_')]) if len(x.rsplit('_')) > 1 else x for x in df['gene_biotype']]
df['gene_name'] = [g.replace('-','.') for g in df['gene_name']]

### gene information
gdf = df[df['feature']=='gene'].copy()
gdf['gene_fullname'] = gdf.apply(lambda x: '_'.join([x['gene_id'],x['gene_name'],x['gene_biotype']]), axis = 1)
fgdf = gdf[['seqname','start','end','strand','gene_fullname']]
fgdf = fgdf.sort_values(by = ['seqname','start'])
fgdf.to_csv(output + '_genes.bed', sep = '\t', index = None)


### finding exons and introns
xdf = {ch: df_ch for ch, df_ch in df.groupby('gene_id')}

def findGeneExonIntron(xdfg):
    xdf_f = {f: df_f for f, df_f in xdfg.groupby('feature')}
    exons_nts = sorted(set([i for idx in xdf_f['exon'].index for i in range(xdf_f['exon'].loc[idx,'start'],xdf_f['exon'].loc[idx,'end'])]))

    i1 = [i for i in range(len(exons_nts)-1) if exons_nts[i+1] != exons_nts[i]+1] + [len(exons_nts)-1]
    i0 = [0] + [i1[i]+1 for i in range(len(i1)-1)]
    if len(i0) > 1 and len(i1) > 1:
        edf = pd.DataFrame([{'start': exons_nts[x0], 'end': exons_nts[x1]+1} for x0, x1 in zip(i0,i1)])
        idf = pd.DataFrame([{'start': exons_nts[x0]+2,'end': exons_nts[x1]-1} for x0, x1 in zip(i1[0:-1],i0[1:])])
        edf['strand'] = xdf_f['gene']['strand'].values[0]
        idf['strand'] = xdf_f['gene']['strand'].values[0]
        edf['gene_fullname'] = '_'.join(xdf_f['gene'][['gene_id','gene_name','gene_biotype']].values[0])
        idf['gene_fullname'] = '_'.join(xdf_f['gene'][['gene_id','gene_name','gene_biotype']].values[0])
        edf['chr'] = xdf_f['gene']['seqname'].values[0]
        idf['chr'] = xdf_f['gene']['seqname'].values[0]
        edf = edf[['chr','start','end','strand','gene_fullname']]
        idf = idf[['chr','start','end','strand','gene_fullname']]
    else:
        edf = pd.DataFrame([{'start': exons_nts[x0], 'end': exons_nts[x1]+1} for x0, x1 in zip(i0,i1)])
        edf['strand'] = xdf_f['gene']['strand'].values[0] 
        edf['chr'] = xdf_f['gene']['seqname'].values[0]
        edf['gene_fullname'] = '_'.join(xdf_f['gene'][['gene_id','gene_name','gene_biotype']].values[0])
        edf = edf[['chr','start','end','strand','gene_fullname']]
        idf = pd.DataFrame(columns = ['chr','start','end','strand','gene_fullname'])
    return (edf, idf)

from multiprocessing import Pool
p = Pool(8)
idf = pd.DataFrame(columns = ['chr','start','end','strand','gene_fullname'])
edf = pd.DataFrame(columns = ['chr','start','end','strand','gene_fullname'])
for mdf in p.imap_unordered(findGeneExonIntron, [xdf[g] for g in xdf.keys()]):
    edf = pd.concat([edf, mdf[0] ])
    idf = pd.concat([idf, mdf[1]])

## clear results
# remove rows where start > end (in the introns i find some)
edf = edf[edf['start']<edf['end']]
idf = idf[idf['start']<idf['end']]

#rep_idf = {i: df_i['gene_fullname'] for i, df_i in idf.groupby(['chr','start','end','strand']) if len(df_i)>1}
#rep_edf = {e: df_e['gene_fullname'] for e, df_e in edf.groupby(['chr','start','end','strand']) if len(df_e)>1}

# pool introns/exons that are identical but have different gene annotations (Gm genes
uidf = [{'chr': i[0], 'start': i[1], 'end': i[2], 'strand': i[3], 'gene_fullname': '-'.join(sorted(np.array(df_i['gene_fullname'])))} for i, df_i in idf.groupby(['chr','start','end','strand'])]
uidf = pd.DataFrame(uidf)[['chr','start','end','strand','gene_fullname']]

uedf = [{'chr': i[0], 'start': i[1], 'end': i[2], 'strand': i[3], 'gene_fullname': '-'.join(sorted(np.array(df_i['gene_fullname'])))} for i, df_i in edf.groupby(['chr','start','end','strand'])]
uedf = pd.DataFrame(uedf)[['chr','start','end','strand','gene_fullname']]

# sort by chrom and start
uedf = uedf.sort_values(by = ['chr','start'])
uidf = uidf.sort_values(by = ['chr','start'])

# print out results
uidf.to_csv(output + '_introns.bed', sep = '\t', index = None,  header=False)
uedf.to_csv(output + '_exons.bed', sep = '\t', index = None, header=False)





