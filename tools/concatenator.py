# Derived from: Mouse gastruloids transcriptomics analysis
# Upstream: https://github.com/anna-alemany/mouseGastruloids_scRNAseq_tomoseq
# Original authors (c) upstream contributors; see upstream history
# Modifications (c) 2025 Moutaz Helal
# License: GPL-3.0
# SPDX-License-Identifier: GPL-3.0-only


# Reads _R1.fq.gz and _R2.fq.gz files, selects reads with proper cell barcode
# and produces a new _cbc.fastq file. Now supports --out to control output path.

import sys, os
import itertools as it
import argparse as argp
import numpy as np
import gzip
import pandas as pd
from pandas.io.parsers import read_csv

#### function to identify cells from barcodes, allowing some edit distances ####
def getCell(bc, bc2cell, hdmax):
    try:
        cell = bc2cell[bc]  # fixed: use the dict passed in
    except:
        cell = 'None'
    return cell

def expandBCset(d, hdmax):
    if hdmax == 0:
        nt = ['N']
        hdmax = 1
    else:
        nt = ['N', 'C', 'T', 'G', 'A']

    dg = {}
    for seq in d:
        try:
            dg[seq].append(d[seq])
        except:
            dg[seq] = [d[seq]]
        i = 0
        while i < hdmax:
            i += 1
            comb = [''.join(l) for l in it.product(nt, repeat = i)]
            for c in comb:
                for p in it.permutations(range(len(seq)), i):
                    s0 = seq
                    for j in range(i):
                        s0 = s0[:p[j]] + c[j] + s0[p[j]+1:]
                    try:
                        dg[s0].append(d[seq])
                    except:
                        dg[s0] = [d[seq]]
    dg2 = dg.copy()
    for seq in dg:
        if len(set(dg[seq])) > 1:
            del dg2[seq]
        else:
            dg2[seq] = dg[seq][0]

    return dg2

#### check input variables ####
parser = argp.ArgumentParser(description = 'Concatenates bcread to bioread qname.')
parser.add_argument('--fqf', help = 'Fastq files names, without _Rx.fastq')
parser.add_argument('--bcread', '-bcr', help = 'read where to find the barcode (umi+cell)', choices = ['R1', 'R2'], default = 'R1')
parser.add_argument('--bioread', '-bior', help = 'read where to find biological information', choices = ['R1', 'R2'], default = 'R2')
parser.add_argument('--lencbc', '-lcbc', help = 'cell barcode length (integer)', type = int, default = 8)
parser.add_argument('--lenumi', '-lumi', help = 'umi length (integer)', type = int, default = 6)
parser.add_argument('--umifirst', help = 'logical variable: umi before cel barcode', action = 'store_true')
parser.add_argument('--cbcfile', '-cbcf', help = 'cell specific barcode file. Please, provide full name')
parser.add_argument('--cbchd', help = 'collapse cell barcodes with the given hamming distance', type = int, default = 0)

# NEW: choose output path prefix (directory + basename, no extension). Script appends _cbc.fastq and .log
parser.add_argument(
    "--out",
    help="Full path for output FASTQ (without extension). Script will append _cbc.fastq",
    default=None
)

args = parser.parse_args()

fqr = args.fqf
bcread = args.bcread
bioread = args.bioread
lcbc = args.lencbc
lumi = args.lenumi
umifirst = args.umifirst
cbcfile = args.cbcfile
hd = args.cbchd

# Determine output prefix (and ensure directory exists)
out_prefix = args.out if args.out else fqr
out_dir = os.path.dirname(out_prefix)
if out_dir:
    os.makedirs(out_dir, exist_ok=True)

#### Define input fastq files ####
fq1 = fqr + '_R1.fq.gz'
fq2 = fqr + '_R2.fq.gz'

if not os.path.isfile(fq1) or not os.path.isfile(fq2):
    print('fastq files not found')
    sys.exit()

#### Read barcodes ####
dbc = read_csv(cbcfile, sep = '\t', index_col=0, header = None)
if not all([len(idx)==lcbc for idx in dbc.index]):
    sys.exit('barcode length provided does not match reference set')
d = {idx: dbc.loc[idx,1] for idx in dbc.index}
bc2sample = expandBCset(d, hd)

#### Do the job ####
out_fastq_path = out_prefix + '_cbc.fastq'
fout = open(out_fastq_path, 'w')

nt = 0
ns = 0
with gzip.open(fq1) as f1, gzip.open(fq2) as f2:
    for idx, (l1, l2) in enumerate(zip(f1, f2)):
        l1, l2 = l1.rstrip().rsplit()[0], l2.rstrip().rsplit()[0]
        l1 = str(l1, 'utf-8')
        l2 = str(l2, 'utf-8')
        l = np.mod(idx,4)
        if l == 0:
            n1, n2 = l1, l2
            if not n1 == n2:
                print (n1, n2)
                sys.exit('fastq files not syncrhonized (@name)')
        if l == 1:
            s1, s2 = l1, l2
        if l == 2:
            p1, p2 = l1[0], l2[0]
            if not p1 == p2: # == '+':
                print(l1, l2)
                print(p1, p2)
                sys.exit('fastq files not synchronized (+)')
        if l == 3:
            q1, q2 = l1, l2
            if len(q1) != len(s1) or len(q2) != len(s2):
                sys.exit('phred and read length not mathch!')

            if bcread == 'R1':
                bcseq = s1[:lumi+lcbc]
                s1 = s1[lumi+lcbc:]
                q1 = q1[lumi+lcbc:]
            elif bcread == 'R2':
                bcseq = s2[:lumi+lcbc]
                s2 = s2[lumi+lcbc:]
                q2 = q2[lumi+lcbc:]
            if not umifirst:
                celbc = bcseq[:lcbc]
                umi = bcseq[lcbc:]
            else:
                celbc = bcseq[lumi:]
                umi = bcseq[:lumi]
            cell = getCell(celbc, bc2sample, hd)
            if cell == 'None':
                continue
            ns += 1
            if bioread == 'R1':
                fout.write('\n'.join([':'.join([n1, bcseq, umi, celbc, str(cell).zfill(3)]), s1, p1, q1, '']))
            elif bioread == 'R2':
                fout.write('\n'.join([':'.join([n2, bcseq, umi, celbc, str(cell).zfill(3)]), s2, p2, q2, ''] ))

nt = (idx+1)/4
fout.close()

#### LOG ####
log_path = out_prefix + '.log'
fout = open(log_path, 'w')
fout.write('=> to generate cbc file <=\n')
fout.write(', '.join(['fastq file:', str(fqr),'\n']))
fout.write(', '.join(['full barcode in:', str(bcread),'\n']))
fout.write(', '.join(['biological read in:', str(bioread), '\n']))
fout.write(', '.join(['cell specific barcode length:', str(lcbc), '\n']))
fout.write(', '.join(['umi length:', str(lumi), '\n']))
fout.write(', '.join(['umi goes first:', str(umifirst),'\n']))
fout.write(', '.join(['total sequenced reads:', str(nt), '\n']))
fout.write(', '.join(['reads with proper barcodes:', str(ns), str(1.0*ns/nt), '\n']))
fout.close()
