#!/bin/sh
for read in *.L1.fastq.gz
  do
    echo $read
    l2=$(echo $read| sed 's/L1/L2/')
    l3=$(echo $read| sed 's/L1/L3/')
    l4=$(echo $read| sed 's/L1/L4/')
    out=$(echo $read| sed 's/.L1//')
    # echo $l2 $l3 $l4 $out
    if [ ! -e merged/$out ]; then
      cat $read $l2 $l3 $l4 > merged/$out
    else
      echo "File exists."
    fi
  done
