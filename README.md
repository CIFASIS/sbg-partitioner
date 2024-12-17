# SBG Partitioner

The goal of this project is to implement a load balancing algorithm to simulate discrete event systems in parallel using set based graphs.


## How to compile and run

`$ make all` or just `$make` compiles and outputs the `sbg-partitioner` binary. To run, it takes these arguments:

* `-f` path to the input file, a json file that represents the model we want to partitionate.
* `-p` number of partitions.
* `-g` [optional argument] output file path.
* `-e` [optional argument] imbalance epsilon, a value between 0 and 1.

`$make sbg-partitioner-metrics` compiles and outputs sbg-partitioner-metrics binary, which partitionates and outputs edge cut, communication volume, maximum communication volume and maximum imbalance of the input and of the files in the directory passed as an argument. To run, it takes these arguments:

* `-f` path to the input file, a json file that represents the model we want to partitionate.
* `-d` path to a directory with txt files that indicate the model partitioning using other algorithms.
* `-p` number of partitions.
* `-e` [optional argument] imbalance epsilon, a value between 0 and 1.

Output files with the metrics will be output in the directory passed as an argument.


`$make sbg-partitioner-exec-time` compiles and outputs sbg-partitioner-exec-time which partitionate the input model 5 times and outputs the average execution time.

* `-f` path to the input file, a json file that represents the model we want to partitionate.
* `-p` number of partitions.
* `-e` [optional argument] imbalance epsilon, a value between 0 and 1.

If the input file is `path/to/file.json`, average execution time will be writen in `path/to/file_${number_of_partitions}_time_exec.txt`.