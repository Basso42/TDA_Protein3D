import re
import statistics
from random import randint
from time import sleep
import statistics
import polars as pl
import pandas as pd
from Bio.Seq import Seq
import requests
import urllib3

import os

import sys
sys.path.append('../')
from src.ESM_functions import fetching_and_writing_3D


if __name__ == '__main__':

    print('Do you want to reset? Y/n')
    choice = str(input())

    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
    
    metadata_path = "/home/onyxia/work/TDA_Protein3D/data/metadata.parquet"
    data_path = "/home/onyxia/work/TDA_Protein3D/data/Tara_relevant_genes.parquet"
    folder_path_3D = "/home/onyxia/work/TDA_Protein3D/data/prot_struct/"


    print('Reading the data')
    relev_genes = pl.read_parquet(data_path)

    #select small proteins and convert their nucleotides in AA
    relev_genes = relev_genes.filter(pl.col('Prot')== True).select(pl.col(['gene', 'sequence']))
    seq2id = relev_genes.with_columns(pl.col('sequence').map_elements(lambda x: str(Seq(x).translate(stop_symbol=''))))
    seq2id = seq2id.with_columns(pl.col("sequence").str.len_bytes().alias("len_prot"))
    small_prots = seq2id.filter(pl.col('len_prot')<=400)
    small_prots = small_prots.drop('len_prot')
    
    small_prots = small_prots.sort(pl.col('gene'))
    print('Cleaning done')
        

    if choice == 'Y':

        #reset folder
        files = glob.glob(struct_path)
        for f in files:
            os.remove(f)

        #reset metadata
        metadata = pd.DataFrame({'Tara_ID': pd.Series(dtype='str'),
                       'Avg_pLDDT': pd.Series(dtype='float')})
        metadata = pl.DataFrame(metadata)
        metadata.write_parquet(metadata_path)

    metada = pl.read_parquet(metadata_path)
    length = small_prots.shape[0] - metada.shape[0]
    missing_3D = small_prots.tail(length+1)

    print(missing_3D.shape)
    
    for i, row in enumerate(missing_3D.rows(named=True)):
        
        seq = row['sequence']
        id = row['gene']
        fetching_and_writing_3D(seq, id, folder_path_3D, metadata_path)
        
        print(i, id)
        