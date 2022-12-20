# This script launches all calculations required for the proof of the
# middle chunk of the interval. It launches several Slurm tasks to
# process different parts of the interval.
# The numbers given here are currently estimates. They should be
# updated once a full run has been done.

# 1:1
# Takes around 44 hours in total. To be conservative the time is set
# to 48 hours.
sbatch --time=48:00:00 -p shared -o PDC/logs/run_proof_1.o -e PDC/logs/run_proof_1.e PDC/scripts/run_proof.sh 1 1

# 2:2
# Takes around 8 hours in total. To be conservative the time is set
# to 12 hours.
sbatch --time=12:00:00 -p shared -o PDC/logs/run_proof_2.o -e PDC/logs/run_proof_2.e PDC/scripts/run_proof.sh 2 2

# 3:9
# Takes around 13 hours in total. To be conservative the time is set
# to 18 hours.
sbatch --time=18:00:00 -p shared -o PDC/logs/run_proof_3.o -e PDC/logs/run_proof_3.e PDC/scripts/run_proof.sh 3 9

# 10:end
# Takes around 16 hours in total. To be conservative the time is set
# to 24 hours.
sbatch --time=24:00:00 -p shared -o PDC/logs/run_proof_4.o -e PDC/logs/run_proof_4.e PDC/scripts/run_proof.sh 10 30
