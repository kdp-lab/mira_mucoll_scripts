Files to run simulation, digitization, and reco are in https://github.com/mlittmann/llp_sim_code. This contains everything for running with MuColl_v1 geometry as well as MAIA, and also (hopefully) contains what one would need for editing Marlin track processors and running locally with those changes.

In this repo are the scripts I've been using for running any sort of analysis. efficiency/efficiency_from_slcio.py runs overall tracking plots, there is also an older version that does this from .jsons but I have found it's most efficient to just run directly with slcios. I've been saving things to caches as well in order to make plots in the same file as the data collection without taking too much time.


