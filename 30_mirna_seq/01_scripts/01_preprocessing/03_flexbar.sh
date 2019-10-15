#!/bin/sh
for read in *.fastq.gz
  do
    echo $read
    out=$(echo $read| sed 's/.fastq.gz//')
    # info=$(echo $read| sed 's/_R1_merge.fq/_info.txt/')
    if [ ! -e flexbar_q/$read ]; then
      flexbar -r $read -a adapters.fa -q TAIL -qf sanger -qw 4 -min-read-length 16 -n 1 --zip-output GZ -t flexbar_q/$out
    else
      echo "File exists."
    fi
  done
