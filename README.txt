Simulation code for "Acceleration by Strong Interactions", https://arxiv.org/abs/1706.07850

List of dependencies:
- cmake version 3.3.0
- gcc (GCC) 4.4.7 20120313 (Red Hat 4.4.7-17)
- numpy 1.11.0
- pxl 3.5: https://forge.physik.rwth-aachen.de/public/pxl/3.1/
- pyROOT 5.34.34



Installation:
1. cd into the directory
2. mkdir build
3. cd build
4. cmake ..
5. make
This installs the program main in that build directory.


Running the trajectory computation from the main directoy:
- cd into the directory
- ./build/trajectory sim_results.csv 100  (simulates 100 trajectories assuming the repulsive potential)

Plotting the angular distances for the attractive and repulsive case from the main directory:
- python relPtPlotter.py sim_results.csv relPt_gain.pdf 


Todo !:
- computation and usage of correct event weights
