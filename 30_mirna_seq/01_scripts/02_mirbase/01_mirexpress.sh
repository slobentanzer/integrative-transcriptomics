#!/bin/sh
mkdir -p mirexpress/results
for read in *.fastq
  do
    echo $read
    dir=$(echo $read| sed 's/.fastq//')
    if [ ! -d mirexpress/results/$dir ]; then
      echo "Create output folder $dir."
      mkdir mirexpress/results/$dir

      echo "Parsing Raw data."
      ~/Programs/miRExpress/src/Raw_data_parse -i $read -o mirexpress/$read.stp2

      # echo "Trimming Adapter~/Programs/miRExpress."
      # ~/Programs/miRExpress/src/Trim_adapter -i example.stp1 -t 3_adaptor.txt -o example.stp2
      # skip, already trimmed

      echo "Generating statistics of reads."
      ~/Programs/miRExpress/src/statistics_reads -i mirexpress/$read.stp2 -o mirexpress/$read.read_statistics.txt


      echo "Aligning reads to precursor miRNA.\nPlease wait.(You can increase the performance according your number of processors/cores! Allow up to n aligning threads \"-u n\")"
      ~/Programs/miRExpress/src/alignmentSIMD -r ~/Programs/miRExpress/data_v2/hsa_precursor.txt -i mirexpress/$read.stp2 -o mirexpress/results/$dir/ -u 8

      echo "Producing results~/Programs/miRExpress.\n"
      ~/Programs/miRExpress/src/analysis -r ~/Programs/miRExpress/data_v2/hsa_precursor.txt -m ~/Programs/miRExpress/data_v2/hsa_miRNA.txt -d mirexpress/results/$dir/ -o read_align_premiRNA.txt -t miRNA_expression.txt -s ~/Programs/miRExpress/data_v2/AllPrecursors_Structures.txt

      echo "*********************\nResults is in folder \"mirexpress/results/$dir/\"\nMicroRNA Expression file is at mirexpress/results/$dir/miRNA_expression.txt\nReads aligned to precursor miRNA files is at mirexpress/results/$dir/read_align_premiRNA.txt\n*********************\n";
    else
      echo "File exists."
    fi
  done
