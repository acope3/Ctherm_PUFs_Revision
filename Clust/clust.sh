#!/bin/bash

~/.local/bin/clust ../WGCNA/Cleaned_data/hpt_cleaned_imputed_data_transposed.tsv -n 4 -o hpt_no_norm/ -r hpt_wgcna_cleaned.txt 
~/.local/bin/clust ../WGCNA/Cleaned_data/ll1210_cleaned_imputed_data_transposed.tsv -n 4 -o LL1210_no_norm/ -r ll1210_wgcna_cleaned.txt 