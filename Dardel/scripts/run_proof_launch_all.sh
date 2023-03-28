# This script launches all calculations required for the proof of the
# middle chunk of the interval. It launches several Slurm tasks to
# process different parts of the interval.

# The first three intervals have a higher memory usage and for that
# reason we use fewer threads to avoid using too much memory. This is
# slower than if we would have used all threads. However it is not
# that much slower due to SMT, only about 20-30% slower it seems.

# 1:1
# Takes around 12 hours in total.
HCW_THREADS=2 sbatch --time=14:00:00 -p main -o PDC/logs/run_proof_1.o -e PDC/logs/run_proof_1.e PDC/scripts/run_proof.sh 1 1

# 2:2
# Takes around 9 hours in total.
HCW_THREADS=2 sbatch --time=11:00:00 -p main -o PDC/logs/run_proof_2.o -e PDC/logs/run_proof_2.e PDC/scripts/run_proof.sh 2 2

# 3:3
# Takes around 6 hours in total.
HCW_THREADS=2 sbatch --time=8:00:00 -p main -o PDC/logs/run_proof_3.o -e PDC/logs/run_proof_3.e PDC/scripts/run_proof.sh 3 3

# 4:10
# Takes around 10 hours in total.
sbatch --time=12:00:00 -p main -o PDC/logs/run_proof_4.o -e PDC/logs/run_proof_4.e PDC/scripts/run_proof.sh 4 10

# 11:14
# Takes around 8 hours in total.
sbatch --time=10:00:00 -p main -o PDC/logs/run_proof_5.o -e PDC/logs/run_proof_5.e PDC/scripts/run_proof.sh 11 14

# 14:end
# Takes around 7 hours in total.
sbatch --time=9:00:00 -p main -o PDC/logs/run_proof_6.o -e PDC/logs/run_proof_6.e PDC/scripts/run_proof.sh 15
