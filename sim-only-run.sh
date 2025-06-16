#!/bin/bash

# Loop through chunks
for i in {0..1}
do
  echo "Running sim analysis for chunk $i..."
  python sim-only.py --chunk $i #--all-events (uncomment if you want to run over multiple events in a chunk)

  # Optional: delay
  sleep 1
done

echo "All chunks have been processed."
