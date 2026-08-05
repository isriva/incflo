#!/bin/bash

# Default parameter file
PARAM_FILE=${1:-params.txt}

# Read parameters from file
while IFS='=' read -r key value; do
    # Skip empty lines and comments
    [[ -z "$key" || "$key" =~ ^#.*$ ]] && continue
    declare "$key=$value"
done < "$PARAM_FILE"

# echo "kmax = $kmax"

# # Run with parameters
# ./main2d.gnu.MPI.ex restart_file=$restart_file kmin=$kmin kmax=$kmax plot_fourier=0  plot_filter=1

# Extract the directory and the checkpoint name
# TARGET_DIR=$(dirname "$restart_file")
# CHK_NAME=$(basename "$restart_file")

# Save the absolute path to your executable so it can be found after cd
EXEC_PATH="$PWD/main2d.gnu.MPI.ex"

# Move into the target directory
cd "$target_dir" || { echo "Failed to cd into $target_dir"; exit 1; }

# Print out the available checkpoint files in this directory
echo "========================================="
echo "Available chk directories in $PWD:"
ls -d chk* 2>/dev/null || echo "  No chk directories found!"
echo "========================================="


# Outer loop: iterate mathematically over the checkpoint numbers
for (( chk_num=chk_start; chk_num<=chk_end; chk_num+=chk_inc )); do
    
    # Reconstruct the checkpoint name (e.g., "chk17000")
    CHK_NAME="chk${chk_num}"
    echo "========================================="
    echo "Processing checkpoint: $CHK_NAME"
    echo "========================================="

    # Inner loop: logarithmic increment using integer arithmetic
    # kmax = (kmax * kmax_inc) / 10
    for (( kmax=kmax_start; kmax<=kmax_end; kmax=(kmax * kmax_inc) / 10 )); do
        # Construct the expected output filename based on your C++ function
        EXPECTED_FILE="filtered_${chk_num}_${kmin}_${kmax}"

        # Check if the file/directory already exists
        if [ -e "$EXPECTED_FILE" ]; then
            echo "  Output $EXPECTED_FILE already exists. Skipping."
        else
            echo "  Running kmin = $kmin and kmax = $kmax on $CHK_NAME"
            
            mpirun -n 8 "$EXEC_PATH" restart_file="$CHK_NAME" kmin="$kmin" kmax="$kmax" plot_fourier=0 plot_filter=1
        fi
        # kmin_high=$(awk "BEGIN {print sqrt($kmax * $kmax + 0.5)}")
        # kmax_high=$((kmax * 100))

        # "$EXEC_PATH" restart_file="$CHK_NAME" kmin="$kmin_high" kmax="$kmax_high" plot_fourier=0 plot_filter=1
    done
done
