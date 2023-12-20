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





#this script requires you to have the gene
def extract_numbers(file_path):
    with open(file_path, 'r') as file:
        content = file.read()

        # Find the pattern "1.00" followed by numbers between 0 and 1
        pattern = re.compile(r'1\.00(?:\s+|\D+)(0(?:\.\d+)?|1(?:\.0+)?)')
        matches = pattern.findall(content)

        # Extract the matched numbers
        numbers = [float(match) for match in matches]

    return numbers
    

def write_pdb_file(pdb_string, output_file):
    with open(output_file + '.pdb', 'w') as f:
        f.write(pdb_string)


def fold_protein_sequence(sequence):
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    data = sequence
    
    headers = {
    'Content-Type': 'application/x-www-form-urlencoded',
}
    response = requests.post(url, headers = headers, data=data, verify=False)
    
    
    return response.text



def fetching_and_writing_3D(prot_seq : str, tara_id : str, folder_path : str, meta_parquet_path : str):
    """
    Args:
        - prot_seq: protein sequence
        - tara_id: tara_ocean_id of the gene

        - folder_path: path where structure is saved
        - meta_parquet_path: path where metadata parquet file stored (pLDDT of the predicted structure linked with TaraID)


    fc: predicts structure of a prot seq using ESMFold API (if length under 400), writes .pdb file, stores confidence prediction in a metadata parquet file
    """

    #predicts structure
    struct = fold_protein_sequence(prot_seq)

    while struct[0:5] == '{"mes':
        sleep(randint(180, 360))
        struct = fold_protein_sequence(prot_seq)
        pdb_path = folder_path + tara_id
        write_pdb_file(struct, pdb_path)

    #saves structure
    pdb_path = folder_path + tara_id
    write_pdb_file(struct, pdb_path)

    print('Structure saved')
    
    #computes average confidence
    indiv_pLDDT = extract_numbers(pdb_path + '.pdb')
    avg_plddt = statistics.fmean(indiv_pLDDT)
             
    #stores metadata
    new_row = pl.DataFrame({'Tara_ID' : tara_id , 'Avg_pLDDT' : avg_plddt})
    metadata = pl.read_parquet(meta_parquet_path)
    metadata = metadata.extend(new_row)
    metadata.write_parquet(meta_parquet_path)
    

    print('Metadata added')