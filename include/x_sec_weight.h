#ifndef X_SEC_WEIGHT_H
#define X_SEC_WEIGHT_H

#include "math.h"
#include "constants.h"
#include "particle.h"
#include <cstdlib>
#include <utility>

std::pair<double, double> compute_phi_theta(double, double, double, double);

void generateKinematic(Particle*, Particle*, Particle*, Particle*);

double eta_q(double);
double a(double);
double kappa(double);
double xpdf(double);


//class to calculate differential cross section for given angles and energy transfers 
//(main code is copied from integral calculation, i.e. synched with Christian)
class xsecWeight{
  public:
    xsecWeight();

    void setSourceMass(double);  //sources are electrons   
    void setTargetMass(double); 
    void setSourceEnergy(double);
    void _updateParams();
    double _x(double, double);
    double _d2sigma_dthetaq_dt(double, double);

  private:
    double m_source;
    double m_target;
    double E_source;

    double s;
    double p;
    double _beta;
    double _gamma;
    double boosted_Ep;
    double boosted_Ee;

};

#endif