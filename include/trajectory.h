#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include <iostream>
#include "event.h"
#include "simResults.h"

std::pair<double, double> runTrajectorySimulation(bool,double,double,double,double,double,double,simResult&);

#endif