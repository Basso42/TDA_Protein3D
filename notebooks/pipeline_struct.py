import re
import statistics


def extract_numbers(file_path):
    with open(file_path, 'r') as file:
        content = file.read()

        # Find the pattern "1.00" followed by numbers between 0 and 1
        pattern = re.compile(r'1\.00(?:\s+|\D+)(0(?:\.\d+)?|1(?:\.0+)?)')
        matches = pattern.findall(content)

        # Extract the matched numbers
        numbers = [float(match) for match in matches]

    return numbers





def fetching_and_writing_3D(prot_seq : string, tara_dict : dict, seq2id: dict, folder_path : string, meta_parquet_path : string):
    """
    Args:
        - prot_seq: protein sequence
        - tara_dict: ordered dictionary linking ordered 

        - folder_path: path where structure is saved
        - meta_parquet_path: path where metadata parquet file stored (pLDDT of the predicted structure linked with TaraID)


    fc: predicts structure of a prot seq using ESMFold API, writes .pdb file, stores confidence prediction in a metadata parquet file
    """

    #predicts structure
    struct = fold_protein_sequence(prot_seq)

    #saves structure
    pdb_path = folder_path + tara_dict[]
    write_pdb_file(struct, pdb_path)

    print('Structure saved')
    
    #computes average confidence
    indiv_pLDDT = extract_numbers(pdb_path)
    avg_plddt = statistics.fmean(indiv_pLDDT)
             
    #stores metadata
    new_row = {'Tara_ID' :  , 'Avg_pLDDT' : avg_plddt}
    metadata = pl.read_parquet(meta_parquet_path)
    metadata = metadata.append(new_row, ignore_index=True)
    metadata.write_parquet(meta_parquet_path)

    print('Metadata added')
    

    