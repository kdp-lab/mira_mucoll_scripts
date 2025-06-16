#!/bin/bash

# Loop through chunks
for i in {0..1}
do
  echo "Running sim analysis for chunk $i..."
  python sim-only.py --chunk $i

  # Optional: delay
  sleep 1
done

echo "All chunks have been processed."
