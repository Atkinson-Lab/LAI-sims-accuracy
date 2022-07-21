# LAI-sims-accuracy
*Code for determining local ancestry inference accuracy*

simulation.sh - general code for simulating admixed individuals with admix-simu https://github.com/williamslab/admix-simu

phasing-jointcall.sh - Filtering sites and phasing of the 1kg-hgdp joint call dataset.

wgs-simulation-rfmix-jointcall.sh - code used for generating all simulated models and RFMix v1 runs, and preparing files for the accuracy calculation.

accuracy.R - code for calculating true positive rates of RFMix calls of simulations, getting counts of miscalls per error mode between ancestry groups, and getting the positions with highest number of miscalls.

accuracy_error_plots.R - code used to generate the accuracy and error modes manuscript figures.

viterbi2msp.R - code to convert RFMix v1 .viterbi output to RFMix v2 .msp output file formats for use with Tractor.

*IDs of reference samples used to generate admixed haplotypes and RFMix Reference:*

amr_simulation.txt - HGDP AMR population IDs used for simulation 

eur_simulation.txt - 1KG IBS population IDs used for simulation

afr_simulation.txt - 1KG YRI population IDs used for simulation


amr_rfmix.txt - HGDP AMR population IDs used as RFMix reference

eur_rfmix.txt - 1KG IBS population IDs used as RFMix reference

afr_rfmix.txt - 1KG YRI population IDs used as RFMix reference

