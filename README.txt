List of dependencies (tested on SLC6, lx3a03):
- cmake version 3.3.0
- gcc (GCC) 4.4.7 20120313 (Red Hat 4.4.7-17)
- numpy 1.11.0
- pyROOT 5.34.34



Installation:
1. cd into the directory
2. mkdir build
3. cd build
4. cmake ..
5. make
This installs the program main in that build directory.



Running the trajectory computation from the main directoy:
1. ./build/main 0 sim_results/no_repulsion.txt 10000  (simulates 10000 trajectories assuming a non repulsive potential, event weights and angular distances of the quarks after e-8 seconds are stored into sim_results/no_repulsion.txt)
2. ./build/main 1 sim_results/with_repulsion.txt 10000  (same as above with repulsive potential)
Note: runtime for 10000 events (after cuts) is about one minute. This corresponds to an enormous speed-up w.r.t. python implementation.



Plotting the angular distances for the attractive and repulsive case from the main directory:
1. python plotting.py sim_results/with_repulsion.txt sim_results/no_repulsion.txt comparison.pdf (order of the files is relevant, produces an histogram for shape comparison between both cases)
