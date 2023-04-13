# Author: Ariana Brenner Clerkin
# Date: 12 Sept. 2022
# conda env: nanopore2

# Import Modules
import os
import numpy as np
import matplotlib.pyplot as plt
import random
import pandas as pd
from setuptools import setup
from setuptools.command.install import install
from IPython.display import Markdown as md
from Bio import SeqIO
import csv
from itertools import product

#Global Variables
fp = '/lustre/fs4/home/amansisido/scratch/RICC-Y/SSRP-Seq-2/fastq_19bp/Attempt1_BLC2fastq/'
path_to_barcodes = '/ru-auth/local/home/abrenner/myscratch/helping_labmates/demultiplexing/'
# Read in file with Names
df = pd.read_csv(path_to_barcodes+'umi_demultiplex_w_names.csv')
barcode_combos = list(df.itertuples(index=False, name=None))

# Read in FastQ and Allocate reads to correct folder 
counter = 0
for read1v2 in range(1,3):# Delete existing fastqs to ensure you do not write into existing files
    print(read1v2, flush = True)
    fastq_to_split = 'Undetermined_S0_R'+str(read1v2)+'_001.fastq'
    data_fp = '/ru-auth/local/home/abrenner/myscratch/helping_labmates/demultiplexing/read'+str(read1v2)+'/'
    for combo in barcode_combos:
        print(data_fp + combo[1], flush = True)
        if not os.path.isfile(data_fp + combo[2]):
            continue
        os.remove(data_fp + combo[2]) 

    # Iterate through fastq of interest
    for seq_record in SeqIO.parse(fp+fastq_to_split, """fastq"""):
        counter = counter + 1
        # Grab the string of interest from the description
        desc_string = seq_record.description[-28:]
        # Iterate through all combos and evaluate in that combo matches the seq_record
        for combo in barcode_combos:
            # Assess start barcode 
            criteria1 = combo[0][:8] == desc_string[:8]
            # Assess end barcode
            criteria2= combo[1]== desc_string[-8:]
            if criteria1 & criteria2:
    #           # If the file does not exist yet, just save the seq_record to a new file with that name
                if not os.path.isfile(data_fp + combo[2]):
                    SeqIO.write(seq_record, data_fp + combo[2], "fastq")
                # If the files does exist, read it in and add onto it
                else:
                    sequences = [i for i in SeqIO.parse(data_fp + combo[2], """fastq""")]
                    sequences.append(seq_record)
                    SeqIO.write(sequences, data_fp + combo[2], "fastq")
    #             outfile.close()
    print(counter)