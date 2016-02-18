SUBDIR = test-data-builder
EXTRA_DIST = coverage.txt gene2refseq_test.txt.bz2 paired-unified.fastq paired1.fastq paired2.fastq python-mapping-detail.txt.bz2 python-mapping-result.txt refseq.rna.bz2 sample.txt sample.txt.bz2 sample.txt.gz sample.txt.lzo sample.txt.xz simulation.fastq.bz2 simulation.txt test-data-builder unique-10gram-not.txt.bz2 unique-10gram.txt.bz2

simulation.fastq: simulation.fastq.bz2
	bzip2 -dc simulation.fastq.bz2 > simulation.fastq
