# Integrative-transcriptomics
This repository contains data and code used in Lobentanzer et al. 2019, "**Integrative transcriptomics reveals sexually dimorphic microRNA control of the cholinergic/neurokine interface in schizophrenia and bipolar disorder**", published in *Cell Reports* on October 15, 2019.
##### In Brief
Lobentanzer et al. show how bioinformatically supported high-throughput techniques such as short RNA sequencing can bridge the gap between traditional molecular interaction studies and purely bioinformatic prediction paradigms in an example focused on disentangling the sexual dimorphism in microRNA regulation of the cholinergic/ neurokine interface in mental disorders.
##### Highlights
+ Single-cell transcriptomes reveal a unique profile of cortical cholinergic neurons
+ Female- and male-derived cells show distinct neurokine-induced miRNA responses
+ Differentially enriched microRNA families constitute a self-organizing network
+ Integrative analysis identifies mir-10/mir-199 regulators of cholinergic function
##### Summary
RNA sequencing analyses are often limited to identifying lowest p value transcripts, which does not address polygenic phenomena. To overcome this limitation, we developed an integrative approach that combines large-scale transcriptomic meta-analysis of patient brain tissues with single-cell sequencing data of CNS neurons, short RNA sequencing of human male- and female-originating cell lines, and connectomics of transcription factor and microRNA interactions with perturbed transcripts. We used this pipeline to analyze cortical transcripts of schizophrenia and bipolar disorder patients. Although these pathologies show massive transcriptional parallels, their clinically well-known sexual dimorphisms remain unexplained. Our method reveals the differences between afflicted men and women and identifies disease-affected pathways of cholinergic transmission and gp130-family neurokine controllers of immune function interlinked by microRNAs. This approach may open additional perspectives for seeking biomarkers and therapeutic targets in other transmitter systems and diseases.
##### Original Article
https://doi.org/10.1016/j.celrep.2019.09.017
##### Organisation
+ 00_ contains Java code used to create and maintain the transcriptional interaction database, *miRNet*
+ 10_ contains R code used in the unbiased meta-analysis of patient brain samples
+ 20_ contains R code used in the single-cell-sequencing-based analysis 
+ 30_ contains R code and bash scripts used in the analysis of small RNA sequencing of LA-N-2 and LA-N-5
+ 40_ contains R code used to define a transcriptional cholinergic system
+ 50_ contains R code used in the subset-specific analyses of patient data
+ 60_ contains R code used to analyse integrative transcriptional interactions of the aforementioned conditions
+ raw_data contains diverse first and third party datasets needed for the analyses
+ working_data is a directory for intermediate data
+ GO.R is a source script for topGO analysis
+ supplement_to_xlsx.R was used to create supplementary CSVs (in spite of the name)
