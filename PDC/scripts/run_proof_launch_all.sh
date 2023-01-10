# This script launches all calculations required for the proof of the
# middle chunk of the interval. It launches several Slurm tasks to
# process different parts of the interval.
# The numbers given here are currently estimates. They should be
# updated once a full run has been done.

# 1:1
# Takes around 12 hours in total. To be conservative the time is set
# to 16 hours.
# It uses more memory than the other versions and for that reason we
# use fewer threads to avoid using too much memory. This is slower
# than if we would have used all threads, but due to SMT it is not
# that much slower.
HCW_THREADS=2 sbatch --time=16:00:00 -p main -o PDC/logs/run_proof_1.o -e PDC/logs/run_proof_1.e PDC/scripts/run_proof.sh 1 1

# 2:2
# Takes around 9 hours in total. To be conservative the time is set
# to 12 hours.
# Similar to the above this uses more memory than the other versions
# and for that reason we use fewer threads to avoid using too much
# memory. This is slower than if we would have used all threads, but
# due to SMT it is not that much slower.
HCW_THREADS=2 sbatch --time=11:00:00 -p main -o PDC/logs/run_proof_2.o -e PDC/logs/run_proof_2.e PDC/scripts/run_proof.sh 2 2

# 3:3
# Takes around 6 hours in total. To be conservative the time is set
# to 8 hours.
# Similar to the above this uses more memory than the other versions
# and for that reason we use fewer threads to avoid using too much
# memory. This is slower than if we would have used all threads, but
# due to SMT it is not that much slower.
HCW_THREADS=2 sbatch --time=8:00:00 -p main -o PDC/logs/run_proof_3.o -e PDC/logs/run_proof_3.e PDC/scripts/run_proof.sh 3 3

# 4:10
# Takes around 10 hours in total. To be conservative the time is set
# to 12 hours.
sbatch --time=12:00:00 -p main -o PDC/logs/run_proof_4.o -e PDC/logs/run_proof_4.e PDC/scripts/run_proof.sh 4 10

# 11:end
# Takes around 16 hours in total. To be conservative the time is set
# to 24 hours.
sbatch --time=24:00:00 -p shared -o PDC/logs/run_proof_5.o -e PDC/logs/run_proof_5.e PDC/scripts/run_proof.sh 10
