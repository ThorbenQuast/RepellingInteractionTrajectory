#ifndef X_SEC_WEIGHT_H
#define X_SEC_WEIGHT_H

#include "math.h"
#include "constants.h"
#include "particle.h"
#include <cstdlib>
#include <utility>


double eta_q(double);
double a(double);
double kappa(double);
double xpdf(double);

double d2sigma_dthetaq_dt(double, double);


#endif