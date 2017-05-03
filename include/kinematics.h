#ifndef KINEMATICS_H
#define KINEMATICS_H

#include "math.h"
#include "constants.h"
#include "kinematics.h"
#include "particle.h"
#include <cstdlib>
#include <utility>


std::pair<double, double> compute_phi_theta(double, double, double, double);
void generateKinematic(Particle*, Particle*, Particle*, Particle*);

#endif