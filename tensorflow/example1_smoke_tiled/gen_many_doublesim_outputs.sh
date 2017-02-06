#!/bin/bash
# genereate 100 sim data sets

# $0 is the script name, $1 id the first ARG, $2 is second...

for (( i=1; i <= 100; i++ ))    
do
    #MANTA_DISABLE_UI=1 /home/hook/thesis_code_final_slim/manta_code/manta/build/manta /home/hook/thesis_code_final_slim/sobel_test_fixed/genDoubleSim.py
    MANTA_DISABLE_UI=1 .../manta ./manta_genSimData.py
done

