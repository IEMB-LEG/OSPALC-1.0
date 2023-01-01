# OSPALC-1.0

Scripts for index selection in https://doi.org/10.1101/2022.10.26.511284.

Usage: bash ospalc-index-selection.sh forward-primer reverse-primer

For example: bash ospalc-index-selection.sh GTGYCAGCMGCCGCGGTAA GGACTACNVGGGTWTCTAAT

This script is used to rank the suitability of some variable components in the long primers, especially the 8nt indices in OSPALC primers. This script is based on short-blast. Depending on blast, selected forward and reverse indices are in forward/reverse.index.txt files, complete primers are in forward/reverse.primer.txt files, index filling rules for different sequencing platform are in forward/reverse.input.txt files.

![image](https://user-images.githubusercontent.com/61352216/210159065-05a74886-b9d5-45c6-96b5-6b2d355d8b0f.png)
