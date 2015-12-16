#! /usr/bin/env python

import sys
import os
import gzip
from glob import glob
import re


def get_tax_lookup(summary_file):
    # Dictionary to lookup taxonomy ID
    _tax_lookup = {}
    lines = (l.strip('\n').split('\t') for l in open(summary_file,'rU'))
    header = lines.next()
    for l in lines:
        asm_id = l[header.index('ftp_path')].split('/')[-1]
        ## asm_id = '%s_%s' % (l[header.index('# assembly_accession')], '_'.join(l[header.index('asm_name')].split()))
        _tax_lookup[asm_id] = l[header.index('taxid')]
    return _tax_lookup

def get_gi_lookup(genbank_file):
    _gi_lookup = {}
    gblines = (l.strip() for l in gzip.open(genbank_file, 'rb') if l.startswith('VERSION'))
    for gbline in gblines:
        x,acc,gi = gbline.split()
        _gi_lookup[acc] = gi.split(':')[1]
    return _gi_lookup


def main(args):
    if args.sumfile == '':
        tax_lookup = get_tax_lookup('%s.assembly_summary.txt' % args.tlevel)
    else:
        tax_lookup = get_tax_lookup(args.sumfile)
    fasta_files = glob('%s/%s/*/latest_assembly_versions/*/*_genomic.fna.gz' % (args.db,args.tlevel))
    for ff in fasta_files:
        print >>sys.stderr, 'Processing %s' % ff
        # Assembly ID
        asm_id = ff.split('/')[-1].split('_genomic')[0]
        # Taxonomy ID
        if asm_id in tax_lookup:
            ti = tax_lookup[asm_id]
        else:
            asm_prefix = asm_id.split('.')[0]
            prematch = [k for k in tax_lookup.keys() if k.startswith(asm_prefix)]
            if len(prematch)==1:
                ti = tax_lookup[prematch[0]]
            else:
                print >>sys.stderr, 'Assembly ID %s is not found, assigning %s' % (asm_id, asm_prefix)
                ti = asm_prefix
        
        # Mapping from accession to gi
        gf = '%s.gbff.gz' % '.'.join(ff.split('.')[:-2])
        assert os.path.exists(gf), "GFF file %s does not exist" % gf 
        gi_lookup = get_gi_lookup(gf)
        
        for l in gzip.open(ff, 'rb'):
            if l.startswith('>'):
                acc = l.split()[0].strip('>')
                gi = gi_lookup[acc]
                print >>args.outfile, '>ti|%s|gi|%s|ref|%s| %s' % (ti, gi, acc, ' '.join(l.split()[1:]))
            else:
                print >>args.outfile, l.strip('\n')

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Get sequences with taxonomy id appended')
    parser.add_argument('--db', help="Database", default="refseq")
    parser.add_argument('--sumfile', help="Summary file", default="")
    parser.add_argument('tlevel', help="Taxonomy level to process")
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout)
    main(parser.parse_args())
