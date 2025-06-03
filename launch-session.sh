#!/bin/bash

# Default values
DEFAULT_TIME="01:00:00"
DEFAULT_CPUS=1
DEFAULT_MEMORY="8G"

usage() {
  echo "Usage: $0 [TIME] [CPUS] [MEMORY] [MOUNT1] [MOUNT2] ..."
  echo "Example: $0 02:00:00 4 16G /home/$USER/data /scratch"
  echo "Defaults:"
  echo "  TIME = $DEFAULT_TIME"
  echo "  CPUS = $DEFAULT_CPUS"
  echo "  MEMORY = $DEFAULT_MEMORY"
}

# If user passes -h or --help, show usage and exit
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
  usage
  exit 0
fi

# Assign args or defaults
TIME=${1:-$DEFAULT_TIME}
CPUS=${2:-$DEFAULT_CPUS}
MEMORY=${3:-$DEFAULT_MEMORY}

# Shift out the first 3 args, now $@ are mounts
shift 3

# Construct mount options for singularity
MOUNTS=""
for dir in "$@"; do
  MOUNTS+=" -B $dir:$dir"
done

echo "Launching SLURM job with:"
echo "  Time: $TIME"
echo "  CPUs: $CPUS"
echo "  Memory: $MEMORY"
echo "  Mounts: $@"

# Run the srun command remotely on hpc7 over ssh:
ssh -t "$(whoami)@hpc7" bash -c "'
  module load singularity
  srun --pty --partition=batch --nodes=1 --ntasks=1 --cpus-per-task=$CPUS --mem=$MEMORY --time=$TIME singularity exec $MOUNTS /home/$(whoami)/CSE_MSE_RXF131/lab-staging/mds3/AdvManu/libslm.sif bash
'"