import json

file_path = "/scratch/miralittmann/analysis/mira_analysis_code/efficiency/bib/loose/3500_10/3500_10_reco0.json"
with open(file_path) as file:
    data = json.load(file)
    print(data.keys())