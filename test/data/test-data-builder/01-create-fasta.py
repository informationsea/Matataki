#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import csv
import sys
import Bio.SeqIO
import commonlib
import os.path
import os

def _main():
    parser = argparse.ArgumentParser(description="Create test data for quickexpression")
    parser.add_argument('-g', '--genes', type=int, nargs="*", default=[12295, 12296, 12297, 12298, 18552, 18553, 18554], help="default: %(default)s")
    parser.add_argument('--output-gene2refseq', type=commonlib.FileType('w'), default=os.path.join(os.getcwd(), os.path.dirname(sys.argv[0]), '../gene2refseq_test.txt.bz2'), help="default: %(default)s")
    parser.add_argument('gene2refseq', type=commonlib.FileType('r'), help="gene2refseq")
    parser.add_argument('refseq_rna', type=commonlib.FileType('r'), nargs="+", help="RefSeq RNA FASTA")
    parser.add_argument('--output-refseq', type=commonlib.FileType('w'), default=os.path.join(os.getcwd(), os.path.dirname(sys.argv[0]), '../refseq.rna.bz2'), help="default: %(default)s")
    options = parser.parse_args()

    available_refseq = set()
    gene2refseq_reader = csv.reader(options.gene2refseq, delimiter='\t', quotechar=None)
    gene2refseq_writer = csv.writer(options.output_gene2refseq, delimiter='\t', quotechar=None)
    gene2refseq_reader.next() # skip header
    for row in gene2refseq_reader:
        if int(row[1]) in options.genes:
            available_refseq.add(row[3])
            gene2refseq_writer.writerow(row)
    options.output_gene2refseq.close()

    for onefile in options.refseq_rna:
        for record in Bio.SeqIO.parse(onefile, 'fasta'):
            if record.id.split('|')[3] in available_refseq:
                Bio.SeqIO.write(record, options.output_refseq, 'fasta')


if __name__ == '__main__':
    _main()
