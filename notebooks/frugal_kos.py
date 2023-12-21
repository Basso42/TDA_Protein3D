import polars as pl
import pandas as pd
import os
from Bio.Seq import Seq
from random import randint
from time import sleep

import sys
sys.path.append('../')
from src.functions import import_s3_parquet
from src.ESM_functions import fetching_and_writing_3D

import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)


def get_layers_for_kos(df_genes_path, list_kos, metadata_samples_path):
    """
    This function takes a parquet file with genes, a list of ko and a path for metadata of Tara samples (table 4) of https://zenodo.org/records/3539258/files/Salazar_et_al_2019_Suppl_Info.xlsx

    and outputs a dataframe containing the corresponding KO_s with an additional column indicating the layer of the sample where the gene comes from
    
    """

    df = pl.read_parquet(df_genes_path) \
           .filter(pl.col("KO").is_in(list_kos)) \ #filter KO
           .with_columns(pl.col('gene').map_elements(lambda x: "_".join(x.split("_", 2)[:2])).alias('sample')) #get id sample
        
    print(df.select(pl.col('sample')))    

    
    layers = ['SRF', 'DCM', 'MES'] #the 3 different layers
    metadata = pl.read_csv(metadata_samples_path)   


    #get all the sample ids for a given layer
    for i, layer in enumerate(layers):
        df_tmp = metadata.filter(pl.col('Layer') == layer)
        
        samples = df_tmp.select(pl.col('PANGAEA sample id')).rows()
        globals()[layer] = [elem[0] for elem in samples] 
       
        
    #add the layer column to metaG df
    df = df.with_columns(pl.when(pl.col("sample").is_in(MES)).then(pl.lit('MES')) \
                            .when(pl.col("sample").is_in(SRF)) \
                            .then(pl.lit('SRF')) \
                            .otherwise(pl.lit('DCM')) \
                            .alias('Layer') 
)


    return df

data_path = "/home/onyxia/work/TDA_Protein3D/data/Tara_relevant_genes.parquet"
path_for_kos = "/home/onyxia/work/TDA_Protein3D/data/3D_for_KOs/"
major_KO_s = ['K19736', 'K00275', 'K05934']




"""
Here we iterate for the list of KO_s we stored in major_KO_s.
For each of them :

    - We get add the layers to the tara dataframe with for each gene
    - We clean the original tara dataset (converting them to proteins and verifying length <400)
    - We predict each gene structure with the ESM API and store it in a specific folder, with the pLDDT transcribed in the metadatafile
    
"""
for ko in major_KO_s:
    print(ko)
    df = get_layers_for_kos(data_path, [ko], "/home/onyxia/work/TDA_Protein3D/data/Table_W4.csv").select(pl.col(["gene","sequence", "Layer"]))
    
    df = df.with_columns(pl.col('sequence').map_elements(lambda x: str(Seq(x).translate(stop_symbol=''))))
    df = df.with_columns(pl.col("sequence").str.len_bytes().alias("len_prot"))


    #check that prot length inferior to 400 (max authorized by ESMFold API)
    small = df.filter(pl.col('len_prot')<=400) 
    print(df.shape[0]- small.shape[0])

    #à changer par la suite (corriger le nom et ranger les fichiers de métadata dans les dossiers correspondants)
    if not os.path.exists(path_for_kos + ko):
        os.makedirs(path_for_kos + ko)
        os.makedirs(path_for_kos + ko + '/struct/')

        #creates metadataparquet
        metadata = pd.DataFrame({'Tara_ID': pd.Series(dtype='str'),
                   'Avg_pLDDT': pd.Series(dtype='float')})
        metadata = pl.DataFrame(metadata)
        metadata.write_parquet(KO_meta_path)

    KO_struct_path = path_for_kos + ko + '/struct/'
    KO_meta_path = path_for_kos + ko + 'metadata.parquet'

    genes_meta = pl.read_parquet(KO_meta_path)
    length = df.shape[0] - genes_meta.shape[0]
    print(length)

    if length > 0: 
        
        df = df.tail(length+1)

        for i, row in enumerate(df.rows(named=True)):
       
           seq = row['sequence']
           id = row['gene']
           fetching_and_writing_3D(seq, id, KO_struct_path, KO_meta_path)
       
           print(i, id)