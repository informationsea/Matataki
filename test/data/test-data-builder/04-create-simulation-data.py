#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import commonlib
import collections
import csv
import os
import os.path
import sys
import random

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

def _main():
    parser = argparse.ArgumentParser(description="Create simulation FASTQ from RefSeq")
    parser.add_argument('--refseq-rna', type=commonlib.FileType('r'), nargs="?", help="default: %(default)s", default=os.path.join(os.getcwd(), os.path.dirname(sys.argv[0]), '../refseq.rna.bz2'))
    parser.add_argument('--unique-ngram', type=commonlib.FileType('r'), nargs="?", help="default: %(default)s", default=os.path.join(os.getcwd(), os.path.dirname(sys.argv[0]), '../unique-10gram.txt.bz2'))
    parser.add_argument('--simulation-fastq', type=commonlib.FileType('w'), nargs="?", help="default: %(default)s", default=os.path.join(os.getcwd(), os.path.dirname(sys.argv[0]), '../simulation.fastq.bz2'))
    parser.add_argument('--simulation-result', type=commonlib.FileType('w'), nargs="?", help="default: %(default)s", default=os.path.join(os.getcwd(), os.path.dirname(sys.argv[0]), '../simulation.txt'))
    parser.add_argument('--fastq-length', type=int, default=50, help="default: %(default)s")
    options = parser.parse_args()

    unique_ngram_reader = csv.reader(options.unique_ngram, delimiter='\t', quotechar=None)
    unique_ngram2gene = dict()
    refseq2gene = dict()
    for row in unique_ngram_reader:
        unique_ngram2gene[row[0]] = row[1]
        for one in row[2].split(', '):
            refseq2gene[one] = row[1]

    ngram = len(unique_ngram2gene.keys()[0])
    if any([len(x) != ngram for x in unique_ngram2gene.iterkeys()]):
        print >>sys.stderr, "Invalid unique ngram file"
        exit(1)

    sequence_names = list()
    sequences = list()
    for record in SeqIO.parse(options.refseq_rna, 'fasta'):
        refseq_id = record.id.split('|')[3]
        sequence_names.append((refseq_id, refseq2gene[refseq_id]))
        sequences.append(str(record.seq))

    simulation_count = [int(2**(random.random()*12)*30) for x in sequences]
    unique_count = [0] * len(sequences)

    fastq_sequences = list()
    
    for c, (oneseq, simcount) in enumerate(zip(sequences, simulation_count)):
        for i in range(simcount):
            fragment_pos = random.randint(0, len(oneseq)-options.fastq_length-1)
            fragment = oneseq[fragment_pos:fragment_pos + options.fastq_length]
            fastq_sequences.append(fragment)
            
            found_gene_id = None
            for j in range(len(fragment) - ngram):
                if fragment[j:j+ngram] in unique_ngram2gene:
                    found_gene_id = unique_ngram2gene[fragment[j:j+ngram]]
            if found_gene_id:
                unique_count[c] += 1
            

    print sequence_names
    print simulation_count
    print unique_count

    random.shuffle(fastq_sequences)
    for i, one in enumerate(fastq_sequences):
        record = SeqRecord(Seq(one, IUPAC.unambiguous_dna),
                           id="SIMULATION"+str(i),
                           description="SIMULATION"+str(i),
                           name="SIMULATION"+str(i))
        record.letter_annotations['phred_quality'] = range(len(one))
        SeqIO.write(record,
                    options.simulation_fastq, 'fastq')

    result_writer = csv.writer(options.simulation_result, delimiter='\t', quotechar=None)
    result_writer.writerow(['Refseq ID', 'GeneID', 'Simulation Count', 'Actual Unique N-Gram Count'])
    for name, simcount, uniqcount in zip(sequence_names, simulation_count, unique_count):
        result_writer.writerow([name[0], name[1], simcount, uniqcount])

if __name__ == '__main__':
    _main()
