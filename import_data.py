import os
import s3fs


data_path = '/home/onyxia/work/TDA_Protein3D/data/'
if not os.path.exists(data_path):
    os.makedirs(data_path)

BUCKET='gamer35/'
folder_s3 = 'KEGG_db/3D_for_KOs'



S3_ENDPOINT_URL = "https://" + os.environ["AWS_S3_ENDPOINT"]
fs = s3fs.S3FileSystem(client_kwargs={'endpoint_url': S3_ENDPOINT_URL})

#get all the predicted structures per ko
fs.get(BUCKET+folder_s3, "/home/onyxia/work/TDA_Protein3D/data/", recursive = True)

#get the partially filtered Tara dataset (has a KO)
fs.get(BUCKET + 'KEGG_db/Tara_relevant_genes.parquet','/home/onyxia/work/TDA_Protein3D/data/')