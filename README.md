now running analysis with condor, so should just have to change reco/sim_analysis.submit. this executes reco/sim_analysis.sh which calls reco/sim_analysis.py. then efficiency.py is used to calculate efficiency
most recent raw .slcio files are saved in /ospool/uc-shared/project/futurecolliders/miralittmann/< reco or sim or digi >/efficiency/< bib or nobib >/< loose or medium or tight >/ < sample >/
and most recent reco .json files are in /scratch/miralittmann/analysis/mira_analysis_code/efficiency/< bib or nobib > /< loose or medium or tight >/< sample >/
