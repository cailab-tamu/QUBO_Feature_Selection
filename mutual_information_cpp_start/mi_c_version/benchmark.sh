#!/bin/bash

# Define the number of threads to test
threads=(1 2 4)

# Define the output file

# Create the output file (overwrite existing file)

# Loop through each thread count
for ((i=0; i<${#threads[@]}; i++)); do
  if [[ ${threads[$i]} -eq 1 ]]; then
    # Run the serial version and capture output
    time ./mi_serial

  else
    # Run the OpenMP version and capture output
   time ./mi_openmp ${threads[$i]} 
  fi

done

