# LAI-sims-accuracy
*Code for determining local ancestry inference accuracy*

simulation.sh - code for simulating admixed individuals with admix-simu https://github.com/williamslab/admix-simu

run_rfmix.sh - code used for RFMix v1 runs of simulations.

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

