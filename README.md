Matataki
========

How to build
------------

1. Run `./autogen.sh` if you cloned from GitHub (This step is not
   required if you downloaded source code from release page.)
2. Run `./configure --disable-shared`
3. Run `make`
4. Run `sudo make install` if you want to install globally
   * Copy `src/matataki`, `src/matataki-builddb`,
     `src/refseqextract` to your directory that is contained in PATH
     variable.

How to use
----------

### Build index file

WARNING: building index requires around 30 GB RAM. We recommend to use 
pre-built index files to quantify human and mouse 

1. Download RefSeq FASTA formatted RNA Sequences from
   [NCBI](ftp://ftp.ncbi.nlm.nih.gov/refseq). (For example,
   ftp://ftp.ncbi.nlm.nih.gov/refseq/M_musculus/mRNA_Prot/mouse.1.rna.fna.gz)
2. Download `refseq2gene` from [NCBI](ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz).
3. Extract and clean sequences from RefSeq.
   * Run `refseqextract -o OUTPUT_FASTA -t TAXID GENE2REFSEQ REFSEQ_FASTA...`
   * For example, run `refseqextract -o clean-mouse.fna.gz -t 10090 gene2refseq.gz *.rna.fna.gz`
     to create a mouse RNA sequence FASTA.
4. Run `matataki-builddb -m -n 32 GENE2REFSEQ REFSEQ_FASTA DBFILE`
   to build index file.
   * For example, run `matataki-builddb -m -n 32 gene2refseq.gz clean-mouse.fna.gz mouse-32.idx`
     to build a mouse index file.

This step is required only once for each species.

### Quantify gene expression

1. Prepare RNA-Seq FASTQ file.
2. Run `matataki -o OUTPUT INDEXDB FASTQ [FASTQ]`
  * For example, run `matataki -o result.txt mouse-32.idx SEQUENCES.fastq`

### Note

* Matataki can read and write gzip or bzip2 compressed file transparently.

License
-------

GNU General Public License version 3 or later

Copyright
---------

Copyright (C) 2016 OKAMURA, Yasunobu
