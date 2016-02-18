#!/bin/bash

cd ${0%/*}

function run() {
  if ! "$@"; then
    echo "ERROR on running" "$@"
    exit 1
  fi
}

DATE=`date +%Y%m%d`

OUTPUT_DIR="download-$DATE"
mkdir -p ${OUTPUT_DIR}

pushd $OUTPUT_DIR

mkdir rawfile

pushd rawfile
wget ftp://anonymous@ftp.hgc.jp/pub/mirror/ncbi/gene/DATA/gene2refseq.gz
#wget 'ftp://ftp.hgc.jp/pub/mirror/ncbi/refseq/H_sapiens/mRNA_Prot/human.*.rna.fna.gz'
wget 'ftp://ftp.hgc.jp/pub/mirror/ncbi/refseq/M_musculus/mRNA_Prot/mouse.*.rna.fna.gz'
#wget 'ftp://ftp.hgc.jp/pub/mirror/ncbi/refseq/release/plant/plant.*.rna.fna.gz'
#wget 'ftp://ftp.hgc.jp/pub/mirror/ncbi/refseq/release/invertebrate/invertebrate.*.rna.fna.gz'
#wget 'ftp://ftp.hgc.jp/pub/mirror/ncbi/refseq/release/vertebrate_other/vertebrate_other.*.rna.fna.gz'
#wget 'ftp://ftp.hgc.jp/pub/mirror/ncbi/refseq/release/fungi/fungi.*.rna.fna.gz'
popd

popd


if [ -e download-current ];then
    rm download-current
fi

ln -s $OUTPUT_DIR download-current
