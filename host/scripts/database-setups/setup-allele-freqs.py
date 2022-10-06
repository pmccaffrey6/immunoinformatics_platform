import requests
import pandas as pd
import numpy as np
import glob
import os
import subprocess
from subprocess import CalledProcessError

FORCE = True
BYPASS = False

allele_database_map = pd.read_csv("../../host-data/allele_frequency_database_mapping.csv")

for idx,row in allele_database_map.iterrows():

    print(row["id"])

    dataset_id = row["id"]
    dataset_name = row["name"]

    if (not BYPASS) and ( (FORCE) or (not os.path.exists(f"../../host-data/allele-frequency-raw-tables/allele_frequency_{dataset_id}.tsv")) ):
        url = f"http://www.allelefrequencies.net/tools/getrawdata.asp?pop_id={dataset_id}&resolved=true"
        try:
            print(dataset_id)
            curl_output = subprocess.check_output(["curl", url])

            curl_output = [i.strip('\r').split(',') for i in curl_output.decode("utf-8").split('\n') if not i in ['\r','']]
            df = pd.DataFrame(curl_output)
            colnames = ["subject_id"] + [i.split("*")[0]+"_"+str(idx) for idx,i in enumerate(df.iloc[0].values[1:])]
            allele_types = sorted(list(set([i.split("_")[0] for i in colnames[1:]])))
            df.columns = colnames

            melts = []
            for allele_type in allele_types:
                melted = pd.melt(df, id_vars=['subject_id'], value_vars=[i for i in colnames if i.startswith(allele_type)])
                melts.append(melted)

            regroup = pd.concat(melts).rename(columns={"variable":"allele_type", "value":"locus"})
            regroup["allele_type"] = regroup["allele_type"].apply(lambda x: x.split("_")[0])
            regroup["dataset_id"] = dataset_id
            regroup["dataset_name"] = dataset_name
            regroup["subject_id"] = [i+"_"+str(dataset_id) for i in regroup["subject_id"]]

            regroup.to_csv(f"../../host-data/allele-frequency-raw-tables/allele_frequency_{dataset_id}.tsv", sep='\t')

            all_tables.append(regroup)

        except:
            pass

if not os.path.exists("../../host-data/allele-frequency-processed-tables"):
    os.makedirs("../../host-data/allele-frequency-processed-tables")

all_tables = [pd.read_csv(df, sep='\t') for df in glob.glob("../../host-data/allele-frequency-raw-tables/allele_frequency_*.tsv")]
pd.concat(all_tables).replace("untyped", np.nan).to_csv(f"../../host-data/allele-frequency-processed-tables/all_allele_frequencies.tsv", sep='\t', index=False)
