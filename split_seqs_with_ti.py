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

def get_taxid_from_genbank(genbank_file):
    ''' Just get the first taxid in genbank file '''
    taxidRE = re.compile('^/db_xref="taxon:(\d+)"')
    gblines = (l.strip() for l in gzip.open(genbank_file, 'rb'))
    for l in gblines:
        m = taxidRE.search(l)
        if m:
            return m.group(1)
    return None

#     
# for rec in split_gb_records(genbank_file):
#         cur_taxid = None
#         cur_acc   = None
#         for l in rec:
#             if l.startswith('VERSION'):
#                 x,acc,gi = l.split()
#                 cur_acc = acc
#             else:
#                 m = taxidRE.search(l)
#                 if m:
#                     cur_taxid = m.group(1)
#                     break
#         print '%s\t%s' % (cur_acc, cur_taxid)
# 
# 

def main(args):
    if args.sumfile == '':
        tax_lookup = get_tax_lookup('%s.assembly_summary.txt' % args.tlevel)
    else:
        tax_lookup = get_tax_lookup(args.sumfile)

    file_index = 1    
    outh = open('%s.%02d.fna' % (args.prefix, file_index), 'w')
    file_char = 0    

    fasta_files = glob('%s/%s/*/latest_assembly_versions/*/*_genomic.fna.gz' % (args.db,args.tlevel))
    for ff in fasta_files:
        # Assembly ID
        asm_id = ff.split('/')[-1].split('_genomic')[0]
        ### Find the taxonomy ID
        if asm_id in tax_lookup:
            # Assembly ID has an exact match in the assembly summary file 
            ti = tax_lookup[asm_id]
        else:
            # Check if Assembly ID prefix matches to assembly summary file
            asm_prefix = asm_id.split('.')[0]
            # print >>sys.stderr, 'Using alternate lookup for %s' % asm_prefix
            prematch = [k for k in tax_lookup.keys() if k.startswith(asm_prefix)]
            if len(prematch)==1:
                ti = tax_lookup[prematch[0]]
            else:
                # Check in genbank file for taxonomy ID
                # print >>sys.stderr, 'Looking for taxid in genbank for %s'  % asm_prefix
                gb_fn = '%s.gbff.gz' % '.'.join(ff.split('.')[:-2])
                assert os.path.exists(gb_fn), "GFF file %s does not exist" % gb_fn
                gb_val = get_taxid_from_genbank(gb_fn)
                if gb_val is not None:
                    ti = gb_val
                    # print >>sys.stderr, 'Found taxid %s for %s' % (gb_val, asm_id)
                else:
                    # Give up, set taxon ID to assembly prefix
                    print >>sys.stderr, 'Assembly ID %s is not found, assigning %s' % (asm_id, asm_prefix)
                    ti = asm_prefix
        
        # Mapping from accession to gi
        gf = '%s.gbff.gz' % '.'.join(ff.split('.')[:-2])
        assert os.path.exists(gf), "GFF file %s does not exist" % gf 
        gi_lookup = get_gi_lookup(gf)

        buffer = []
        buffer_char = 0
        for l in gzip.open(ff, 'rb'):
            if l.startswith('>'):
                acc = l.split()[0].strip('>')
                gi = gi_lookup[acc]
                # print >>args.outfile, '>ti|%s|gi|%s|ref|%s| %s' % (ti, gi, acc, ' '.join(l.split()[1:]))
                headerline = '>ti|%s|gi|%s|ref|%s| %s' % (ti, gi, acc, ' '.join(l.split()[1:]))
                if (file_char + buffer_char) > args.maxchar:
                    print >>sys.stderr, '%s.%02d.fna has %d characters' % (args.prefix, file_index, file_char)
                    # create a new file                
                    outh.close()
                    file_index += 1
                    outh = open('%s.%02d.fna' % (args.prefix, file_index), 'w')
                    file_char = 0
                
                # Write buffer to file
                if len(buffer):                
                    print >>outh, '\n'.join(buffer)
                    file_char += buffer_char
                # Create new buffer
                buffer = [ headerline ]
                buffer_char = 0
            else:
                #print >>args.outfile, l.strip('\n')
                buffer.append(l.strip('\n'))
                buffer_char += (len(l) - 1)

        # Clear buffer for file
        if (file_char + buffer_char) > args.maxchar:
            print >>sys.stderr, '%s.%02d.fna has %d characters' % (args.prefix, file_index, file_char)        
            # create a new file                
            outh.close()
            file_index += 1
            outh = open('%s.%02d.fna' % (args.prefix, file_index), 'w')
            file_char = 0            
        # Write buffer to file
        print >>outh, '\n'.join(buffer)
        file_char += buffer_char
    
    print >>sys.stderr, '%s.%02d.fna has %d characters' % (args.prefix, file_index, file_char)  
    outh.close()

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Get sequences with taxonomy id appended')
    parser.add_argument('--db', help="Database", default="refseq")
    parser.add_argument('--sumfile', help="Summary file", default="")
    parser.add_argument('--maxchar', type=float,  help="Max characters in file", default=4e9)
    parser.add_argument('tlevel', help="Taxonomy level to process")
    parser.add_argument('prefix', help="Output prefix")
    main(parser.parse_args())
