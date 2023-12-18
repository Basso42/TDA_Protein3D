import s3fs
import os
import glob



def import_s3_parquet(s3_folder_path, local_path):
    """
    Imports all the .parquet file to {local_path} from a given S3 folder"
    """
    S3_ENDPOINT_URL = "https://" + os.environ["AWS_S3_ENDPOINT"]
    fs = s3fs.S3FileSystem(client_kwargs={'endpoint_url': S3_ENDPOINT_URL})
    BUCKET='gamer35/'

    file_list = fs.ls(BUCKET + s3_folder_path)
    
    # Loop through the CSV files and download them
    for file_path in file_list:
        if file_path.endswith('.parquet'):
            # Construct the local file path where the CSV will be saved
            local_file_path = os.path.join(local_path, os.path.basename(file_path))
            
            # Download the CSV file from S3 and save it locally
            with fs.open(file_path, 'rb') as s3_file, open(local_file_path, 'wb') as local_file:
                local_file.write(s3_file.read())
    
    print(f"Downloaded parquet files to {local_path}")