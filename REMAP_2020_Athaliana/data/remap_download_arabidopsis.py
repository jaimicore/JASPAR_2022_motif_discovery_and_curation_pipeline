import requests
import os
import re
import pandas

r = requests.get('http://remap.univ-amu.fr:80/api/v1/datasets/findByTaxid/taxid=3702')
json_data = r.json()

dsets_names_list = list()
dsets_TFs_list = list()

for dataset in json_data["datasets"]:
    os.mkdir(dataset["dataset_name"])
    file = dataset["bed_url"]
    file = re.sub("tair10/", "tair10/tf/", file) # needs to be done because API at the moment returns broken URL
    download = requests.get(file)
    open(os.path.join(dataset["dataset_name"], os.path.basename(dataset["bed_url"])), 'wb').write(download.content)
    dsets_names_list.append(dataset["dataset_name"])
    dsets_TFs_list.append(dataset["target_name"])

mapping = pandas.DataFrame(list(zip(dsets_names_list, dsets_names_list, dsets_TFs_list)), columns = ["Dataset_name", "Dataset_name", "TF"])
mapping.to_csv("mapping.tsv", sep = "\t", index = False, header = False)
