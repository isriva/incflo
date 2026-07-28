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
TARGET_DIR=$(dirname "$restart_file")
CHK_NAME=$(basename "$restart_file")

# Save the absolute path to your executable so it can be found after cd
EXEC_PATH="$PWD/main2d.gnu.MPI.ex"

# Move into the target directory
cd "$TARGET_DIR" || exit

# Run with parameters in a loop
for (( kmax=kmax_start; kmax<=kmax_end; kmax+=kmax_inc )); do
    echo "kmax = $kmax"
    # ./main2d.gnu.MPI.ex restart_file="$restart_file" kmin="$kmin" kmax="$kmax" plot_fourier=0 plot_filter=1
    mpirun -n 8 "$EXEC_PATH" restart_file="$CHK_NAME" kmin="$kmin" kmax="$kmax" plot_fourier=0 plot_filter=1

    # kmin_high=$(awk "BEGIN {print sqrt($kmax * $kmax + 0.5)}")
    # kmax_high=$((kmax * 100))

    # "$EXEC_PATH" restart_file="$CHK_NAME" kmin="$kmin_high" kmax="$kmax_high" plot_fourier=0 plot_filter=1
done