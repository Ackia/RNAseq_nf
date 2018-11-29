#!/usr/bin/env Rscript
repo='https://ftp.acc.umu.se/mirror/CRAN/'
install.packages("RMySQL", repos=repo)
install.packages("BiocManager", repos=repo)
BiocManager::install("GenomicFeatures")
BiocManager::install("DESeq2")
BiocManager::install("tximport")
