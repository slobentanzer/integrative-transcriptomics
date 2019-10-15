#!/bin/sh
find . -name "*_L001_*" -type f -exec rename 's/_L001_R1_001.fastq.gz/.L1.fastq.gz/' {} +
find . -name "*_L002_*" -type f -exec rename 's/_L002_R1_001.fastq.gz/.L2.fastq.gz/' {} +
find . -name "*_L003_*" -type f -exec rename 's/_L003_R1_001.fastq.gz/.L3.fastq.gz/' {} +
find . -name "*_L004_*" -type f -exec rename 's/_L004_R1_001.fastq.gz/.L4.fastq.gz/' {} +
