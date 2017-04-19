#ifndef FORCE_H
#define FORCE_H

#include <vector>
#include "constants.h"
#include "particle.h"
#include "pxl/core.hh"

pxl::Basic3Vector Force(Particle*, Particle*, bool);

#endif