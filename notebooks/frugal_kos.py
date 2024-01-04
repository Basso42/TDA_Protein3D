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

    #filter KO and get id sample
    df = pl.read_parquet(df_genes_path) \
           .filter(pl.col("KO").is_in(list_kos)) \
           .with_columns(pl.col('gene').map_elements(lambda x: "_".join(x.split("_", 2)[:2])).alias('sample')) 
        
    #print(df.select(pl.col('sample')))    

    
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
major_KO_s = ['K19736', 'K00275', 'K05934', 'K09023', 'K06871', 'K00331', 'K22024', 'K03214', 'K12064', 'K21947'] #top 4 des trois SHAP des METAT

#


"""
Here we iterate for the list of KO_s we stored in major_KO_s.
For each of them :

    - We add a column to the tara dataframe containing layer of sample for each gene
    - We clean the original tara dataset (converting them to proteins and verifying length <400) /!\ Not verified proteins ! 
    - We predict each gene structure with the ESM API and store it in a specific folder, with the pLDDT transcribed in the metadatafile
    - If some proteins are too long, we write their tara_gene_id in a .txt file to predict them later with ColabFold
    
"""
if __name__ == '__main__':
    for ko in major_KO_s:
        print(ko)
    
        ko_path = path_for_kos + ko
        KO_struct_path = ko_path + '/struct/'
        KO_meta_path = ko_path + f'/{ko}_metadata.parquet'
        
        df = get_layers_for_kos(data_path, [ko], "/home/onyxia/work/TDA_Protein3D/data/Table_W4.csv").select(pl.col(["gene","sequence", "Layer"]))
        
        df = df.with_columns(pl.col('sequence').map_elements(lambda x: str(Seq(x).translate(stop_symbol=''))))
        df = df.with_columns(pl.col("sequence").str.len_bytes().alias("len_prot"))
    
    
        #check that prot length inferior to 400 (max authorized by ESMFold API)
        small = df.filter(pl.col('len_prot')<=400) 
        print("Long proteins : ", df.shape[0]- small.shape[0])
    
        if not os.path.exists(path_for_kos + ko):
            os.makedirs(path_for_kos + ko)
            os.makedirs(path_for_kos + ko + '/struct/')
    
            #creates metadataparquet
            metadata = pd.DataFrame({'Tara_ID': pd.Series(dtype='str'),
                       'Avg_pLDDT': pd.Series(dtype='float')})
            metadata = pl.DataFrame(metadata)
            metadata.write_parquet(KO_meta_path)
    
        #long proteins to predict with alpha_fold
        if df.shape[0]- small.shape[0]>0:
            long_prots = df.filter(pl.col('len_prot')>400) \
                           .select(pl.col('gene'))
            ids = long_prots.rows()
            ids = [elem[0] for elem in ids]
            
            with open(f'{ko_path}/long_prot.txt', 'w') as fp:
                for item in ids:
                    # write each item on a new line
                    fp.write("%s\n" % item)
                print('Long proteins written')        
                                   
    
    
        genes_meta = pl.read_parquet(KO_meta_path)
        length = small.shape[0] - genes_meta.shape[0]
        print("Number of proteins still to fold: ", length)
    
        if length > 0: 
            
            small = small.tail(length+1)
    
            for i, row in enumerate(small.rows(named=True)):
           
               seq = row['sequence']
               id = row['gene']
               fetching_and_writing_3D(seq, id, KO_struct_path, KO_meta_path)
           
               print(i, id)