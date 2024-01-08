#  TDA_Protein3D

This repository contains some scripts that fetch protein structure predictions from the ESM Fold API (such as /notebooks/frugal_kos.py) for certain amino-acid sequences translated from the Tara Ocean prokaryote dataset. 

- /notebooks/frugal_kos.py : predicts and stores 3D structures for proteins (not 100% sure if proteins) having a certain KO (inform the list of KO before running it)

- /notebooks/pipeline_struct.py : predicts and stores 3D structures for all observations in Tara datasets


It also contains a notebook that tests different Topological Data Analyses (/notebooks/Exploring.ipynb) approach on these structures to find some recurrent patterns and study the topological variance of the proteins between layers. 

```import_data.py``` must be run at the beginning of the session on the SSPLab cluster. It downloads the folder containing protein structures for KO's in the /data folder. It is then structured as the following 

```
├── data
│   ├── 3D_for_KOs
│       ├── KOO275
            ├── struct -> *containing the 3D structures of all the Tara ocean proteins under 400 AA with a .pdb file for each protein.*
            ├── KOO275_metadata.parquet -> *contains the id of the gene and the pLLDT confidence for the predicted structure*
            ├── long_prot.txt -> *might not exist, contains the ids of the long proteins that could not be folded (>400 AA)*
        ├── KOO331
        ├── ...
```
