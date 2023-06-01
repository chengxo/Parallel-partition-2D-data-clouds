# Parallel-partition-2D-data-clouds
To compile partition.c, you should install p4est package and compile with it. See the official site of p4est package: https://p4est.github.io/.
For convenience, the small executable "partition" is already included in the repo.

To run the executable, a dataset file and the total number of points needs to be provided as arguments. See sbatch files as examples.
The dataset file needs to separate each coordinate by '\n'.
