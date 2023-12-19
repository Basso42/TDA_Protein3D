### Downloading the data from Luis google drive (preferably with ethernet) ###

# OM-RGC_v2.tsv https://drive.google.com/file/d/1YUNHJXFJ3T8bcoi1YEI8gZ5e6mBzY93s/view?usp=sharing
# To download gene_profile_metaG.tsv from my inria account


"""
Here, find way to run wget in .py file
"""

!wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1YUNHJXFJ3T8bcoi1YEI8gZ5e6mBzY93s' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1YUNHJXFJ3T8bcoi1YEI8gZ5e6mBzY93s" -O OM-RGC_v2.tsv && rm -rf /tmp/cookies.txt

# request creation
q = (
    pl.scan_csv("/home/onyxia/work/TDA_Protein3D/notebooks/OM-RGC_v2.tsv", separator='\t') # Lecture lazy, put your download path here
    .filter(pl.col("KO").is_in(KO))
    .with_columns(pl.col('sequence').str.starts_with('ATG').alias('Prot'))
)

# request execution
df = q.collect()