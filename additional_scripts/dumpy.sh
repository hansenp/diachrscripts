#!/bin/bash
# Usage examples
# dumpy.sh  MIF_R1 MIF_R1 "ERR436029"
# dumpy.sh  MIF_R2 MIF_R2 "ERR436028;ERR436030;ERR436033"
# dumpy.sh  MIF_R3 MIF_R3 "ERR436031;ERR436026"
OUT_DIR=$1
OUT_PREFIX=$2
SRR_LIST=$3
# Init target results
mkdir $OUT_DIR
> $OUT_DIR/$OUT_PREFIX\_1.fastq.gz # Forward
> $OUT_DIR/$OUT_PREFIX\_2.fastq.gz # Reverse
> $OUT_DIR/$OUT_PREFIX\_md5.txt
# Iterate through list of identifiers
IFS=";"
for SRR in $SRR_LIST;
do
  # Download and append
  SRR_TRUNC=$(echo $SRR | cut -c1-6)
  # Forward
  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$SRR_TRUNC/$SRR/$SRR\_1.fastq.gz -O $OUT_DIR/$SRR\_1.fastq.gz
  md5sum $OUT_DIR/$SRR\_1.fastq.gz >> $OUT_DIR/$OUT_PREFIX\_md5.txt
  cat $OUT_DIR/$SRR\_1.fastq.gz >> $OUT_DIR/$OUT_PREFIX\_1.fastq.gz
  rm $OUT_DIR/$SRR\_1.fastq.gz
  # Reverse
  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$SRR_TRUNC/$SRR/$SRR\_2.fastq.gz -O $OUT_DIR/$SRR\_2.fastq.gz
  md5sum $OUT_DIR/$SRR\_2.fastq.gz >> $OUT_DIR/$OUT_PREFIX\_md5.txt
  cat $OUT_DIR/$SRR\_2.fastq.gz >> $OUT_DIR/$OUT_PREFIX\_2.fastq.gz
  rm $OUT_DIR/$SRR\_2.fastq.gz
done
# Compress concatenated files
md5sum $OUT_DIR/$OUT_PREFIX\_1.fastq.gz >> $OUT_DIR/$OUT_PREFIX\_md5.txt
md5sum $OUT_DIR/$OUT_PREFIX\_2.fastq.gz >> $OUT_DIR/$OUT_PREFIX\_md5.txt
