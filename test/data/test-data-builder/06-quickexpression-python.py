#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import commonlib
import csv
import collections
from Bio import SeqIO

def _main():
    parser = argparse.ArgumentParser(description="Python implementation of quick expression")
    parser.add_argument('index', help='tab delimited file', type=commonlib.FileType('r'), nargs='?', default='../unique-10gram.txt.bz2')
    parser.add_argument('fastq', help='FASTQ file', type=commonlib.FileType('r'), nargs='?', default='../simulation.fastq.bz2')
    parser.add_argument('coverage', help='Coverage file', type=commonlib.FileType('r'), nargs='?', default='../coverage.txt')
    parser.add_argument('-o', '--output', help='output path', type=commonlib.FileType('w'), default='../python-mapping-result.txt')
    parser.add_argument('-d', '--detail', help='detail output path', type=commonlib.FileType('w'), default='../python-mapping-detail.txt.bz2')
    options = parser.parse_args()

    # load database
    indexdb = {}
    ngram = 0
    for row in csv.reader(options.index, delimiter='\t', quotechar=None):
        indexdb[row[0]] = row[1]
        ngram = len(row[0])

    transcript_coverage = collections.defaultdict(list)
    reader = csv.reader(options.coverage, delimiter='\t', quotechar=None)
    reader.next() # skip header
    for row in reader:
        print row[0], row[5]
        transcript_coverage[row[0]].append(int(row[5]))

    foundFragments = {x: float(sum(y))/len(y) for (x, y) in transcript_coverage.iteritems()}
    for k, v in foundFragments.iteritems():
        print "Found Fragments: ", k, v, transcript_coverage[k], len(transcript_coverage[k])

    result = collections.defaultdict(int)

    count_fastq = 0
    for record in SeqIO.parse(options.fastq, 'fastq'):
        seq = str(record.seq)
        mapping = ''
        candidate = None
        for i in range(len(seq) - ngram + 1):
            fragment = seq[i:i+ngram]
            if fragment in indexdb:
                if not candidate:
                    candidate = indexdb[fragment]
                    mapping += 'O'
                elif candidate != indexdb[fragment]:
                    candidate = "!"
                    mapping += 'X'
                else:
                    mapping += 'O'
            else:
                mapping += '.'
                    
        if candidate:
            count_fastq += 1
            result[candidate] += 1
            print >>options.detail, '{}\t{}\t{}'.format(record.name, candidate, mapping)
        else:
            result['UNMAPPED'] += 1
            print >>options.detail, '{}\t{}\t{}'.format('UNMAPPED', candidate, mapping)


    fpkm = {x: float(y)/foundFragments[x]*1000/count_fastq*1000000 for x, y in result.iteritems() if x != 'UNMAPPED'}


    for k, v in result.iteritems():
        print >>options.output, '{0}\t{1}\t{2:.22}'.format(k, v, fpkm.get(k, '0.0'))

if __name__ == '__main__':
    _main()
