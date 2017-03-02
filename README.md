# HiC
This repository contains a collection of tools for analysing HiC data<br><br>
<b>ChrInv.py:</b> Determines an average expected number of read for the interaction between two sites of a given distance and uses this to determine sites that are out of place by means of higher than expected interactions, i.e. chromosome inversions<br>
<b>hicpro2heatmap.py:</b> Take a matrix output file from the HiC-Pro pipeline and generates a heatmap for a given chromosome. Can output obs, exp and o/e heatmaps but does not do matrix balancing on raw data<br>
<b>hicup2juicebox.py:</b> Uses the bam file output from the HiCUP pipeline to generate the required sorted text file required for the Juicebox CLT pre command to make a binary .hic file<br>