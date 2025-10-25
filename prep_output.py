import os
import pandas as pd
import glob


def prep_output(file_list, path):
    merged_df= None
    for file in file_list:

        sample_name = os.path.basename(os.path.join(path,file)).replace("_counts.txt", "")

        df = pd.read(file, sep="\t", comment="#")
        counts= df.columns[-1]
        df = [["Geneid", counts]].rename(columns={counts:sample_name})

        if merged_df is None:
            merged_df = df
        else: 
            merged_df = merged_df.merge(df, on="Geneid", how="outer")
    merged_df.to_csv(os.path.join(path, "counts_all_replicates.csv"), sep=",", index=False)

path = "./Data"
file_list = glob.glob(os.path.join(path, "*_counts.txt"))
prep_output(file_list, path)

                     
        
    