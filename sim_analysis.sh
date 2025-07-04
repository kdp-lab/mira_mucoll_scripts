#!/bin/bash
echo $HOSTNAME
echo "<<<Singularity ENVIRONMENT:" $SINGULARITY_NAME
echo "<<<Setup some environment"
echo "source /opt/setup_mucoll.sh --> "
source /opt/setup_mucoll.sh
echo ">>>completed"
echo "<<<Check if we can find executables"
which ddsim
which Marlin
echo "<<<Check if input files were copied from the origin"
ls -lta 

input_file=""
output_directory=""
chunks=""
n_events=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        --inputFile)
            input_file="$2"
            shift 2
            ;;
        --chunks)
            chunks="$2"
            shift 2
            ;;
        --nEvents)
            n_events="$2"
            shift 2
            ;;
        --pId)
            proc_id="$2"
            shift 2
            ;;
        *)
            usage
            ;;
    esac
done

# Construct the nohup command
command="python sim_analysis.py --chunk ${proc_id} --sample ${input_file}"

# Print the constructed command
echo "Executing command: $command"

# Run the command
eval $command

# Copy finished sim file

echo "<<<Delete input files so they don't get transfered twice on exit"
rm -rf sim_condor_analyzer.py
# rm -rf ${input_file}_sim${proc_id}.slcio
echo ">>> Deletions complete. Test job complete"