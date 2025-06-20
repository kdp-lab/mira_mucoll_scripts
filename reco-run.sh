#!/bin/bash

# Loop through chunks
for i in {0..499}
do
  echo " "
  echo "-----------"
  echo "Running chunks_analyzer.py for chunk $i..."
  python reco_analyzer.py --chunk $i

  # Optional: delay
  sleep 1
done

echo "All chunks have been processed."

#echo "Now running aggregator script..."

#python3 aggregate_stau_counts.py

