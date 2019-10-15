#!/bin/sh
mkdir length_filtered
for read in *.fastq
  do
    echo $read
    # out=$(echo $read| sed 's/.fq/_read_lengths.txt/')
    if [ ! -e length_filtered/$read ]; then
      awk 'BEGIN {FS = "\t" ; OFS = "\n"} {
        header = $0 ;
        getline seq ;
        getline qheader ;
        getline qseq ;
        if (length(seq) >= 16 && length(seq) <= 25) {print header, seq, qheader, qseq}
      }' \
      < $read > length_filtered/$read
    else
      echo "File exists."
    fi
  done
