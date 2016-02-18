#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import commonlib
import collections
import csv
import os
import os.path
import sys

from Bio import SeqIO
from Bio import Seq

BASEDIR = os.path.dirname(__file__)

def _main():
    parser = argparse.ArgumentParser(description="Search Unique Sequence")
    parser.add_argument('--gene2refseq', type=commonlib.FileType('r'), help="default: %(default)s", nargs='?', default=os.path.join(BASEDIR, '../gene2refseq_test.txt.bz2'))
    parser.add_argument('--refseq_rna', type=commonlib.FileType('r'), nargs="?", help="default: %(default)s", default=os.path.join(BASEDIR, '../refseq.rna.bz2'))
    parser.add_argument('--unique-ngram', type=commonlib.FileType('w'), nargs="?", help="default: %(default)s", default=os.path.join(BASEDIR, '../unique-10gram.txt.bz2'))
    parser.add_argument('--not-unique-ngram', type=commonlib.FileType('w'), nargs="?", help="default: %(default)s", default=os.path.join(BASEDIR, '../unique-10gram-not.txt.bz2'))
    parser.add_argument('--fragment-specific-unique-ngram', type=commonlib.FileType('w'), nargs="?", help="default: %(default)s", default=os.path.join(BASEDIR, '../unique-10gram-fragment-specific.txt.bz2'))
    parser.add_argument('--ngram', type=int, default=10)
    parser.add_argument('--coverage-output', type=commonlib.FileType('w'), help="default: %(default)s", default=os.path.join(BASEDIR, '../coverage.txt'))
    options = parser.parse_args()

    # load gene2refseq
    gene2refseq = collections.defaultdict(set)
    refseq2gene = dict()

    for row in csv.reader(options.gene2refseq, delimiter='\t', quotechar=None):
        if row[3] == '-': continue
        
        gene2refseq[row[1]].add(row[3])
        refseq2gene[row[3]] = row[1]

    print gene2refseq

    refseq2sequence = dict()
    ngram2refseq = collections.defaultdict(set)

    available_gene2refseq = collections.defaultdict(list)

    # load sequences and calculate uniqueness

    for record in SeqIO.parse(options.refseq_rna, 'fasta'):
        refseq_id = record.id.split('|')[3]
        refseq2sequence[refseq_id] = {'id': refseq_id, 'seq': str(record.seq), 'coverage': list(['.']*len(record.seq)), 'gene': 0, 'fragment': 0}
        available_gene2refseq[refseq2gene[refseq_id]].append(refseq_id)
        for i in xrange(len(record.seq)-options.ngram+1):
            fragment = record.seq[i:i+options.ngram]
            
            ngram2refseq[str(fragment)].add((refseq_id, refseq2gene[refseq_id]))
            ngram2refseq[str(Seq.reverse_complement(fragment))].add((refseq_id, refseq2gene[refseq_id]))

    # write out gene unique and non-unique, isoform specific n-grams

    unique_ngrams = []
    unique_no_common_ngrams = []
    output = csv.writer(options.unique_ngram, delimiter='\t', quotechar=None)
    output_notunique = csv.writer(options.not_unique_ngram, delimiter='\t', quotechar=None)
    output_specific = csv.writer(options.fragment_specific_unique_ngram, delimiter='\t', quotechar=None)
    for k, v in ngram2refseq.iteritems():
        geneset = set([x[1] for x in v])
        refseqset = set([x[0] for x in v])
        if len(geneset) == 1:
            geneid = geneset.pop()
            print refseqset, available_gene2refseq[geneid]
            if refseqset == set(available_gene2refseq[geneid]):
                unique_ngrams.append((k, geneid, ', '.join(refseqset)))
                output.writerow([k, geneid, ', '.join(refseqset)])
            else:
                unique_no_common_ngrams.append((k, geneid, ', '.join(refseqset)))
                output_specific.writerow([k, geneid, ', '.join(refseqset)])
        else:
            output_notunique.writerow([k, ', '.join(geneset), ', '.join(refseqset)])

    # isoform specific reads
    for oneunique in unique_no_common_ngrams:
        #print oneunique
        ngram = len(oneunique[0])
        for refseq in oneunique[2].split(', '):
            record = refseq2sequence[refseq]
            refseq2sequence[refseq]['gene'] = oneunique[1]
            pos = 0
            while True:
                pos = record['seq'].find(oneunique[0], pos)
                if pos < 0: break
                record['coverage'][pos] = '_'
                for i in xrange(ngram-1):
                    if record['coverage'][pos+1+i] == '.':
                        record['coverage'][pos+1+i] = ','
                pos += 1


    # calculate coverage
    for oneunique in unique_ngrams:
        #print oneunique
        ngram = len(oneunique[0])
        for refseq in oneunique[2].split(', '):
            record = refseq2sequence[refseq]
            refseq2sequence[refseq]['gene'] = oneunique[1]
            refseq2sequence[refseq]['fragment'] += 1
            pos = 0
            while True:
                pos = record['seq'].find(oneunique[0], pos)
                if pos < 0: break
                record['coverage'][pos] = 'X'
                for i in xrange(ngram-1):
                    if record['coverage'][pos+1+i] != 'X':
                        record['coverage'][pos+1+i] = '*'
                pos += 1

    # write out coverage
    coverage_result = csv.writer(options.coverage_output, delimiter='\t', quotechar=None)
    coverage_result.writerow(['GeneID', 'RefSeqID', 'Length', '# of covered base', '# of unique fragment', '# of found fragmnet', 'coverage', 'cover'])
    for k, v in refseq2sequence.iteritems():
        coverage = v['coverage']
        covered = len([x for x in v['coverage'] if x == 'X' or x == '*'])
        starts = len([x for x in v['coverage'] if x == 'X'])
        assert len(coverage) == len(v['seq'])
        percent = starts/float(len(coverage)-ngram+1)
        coverage_result.writerow([v['gene'], k, len(coverage), covered, v['fragment'], starts, str(percent), ''.join(coverage)])

if __name__ == '__main__':
    _main()
