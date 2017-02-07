#!/bin/bash
# genereate 10 sim data sets

for (( i=1; i <= 10; i++ ))    
do
    MANTA_DISABLE_UI=1 ../build/manta ./manta_genSimData.py
done

