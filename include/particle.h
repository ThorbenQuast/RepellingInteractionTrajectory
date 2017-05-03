#ifndef PARTICLE_H
#define PARTICLE_H


#include "constants.h"
#include <vector>
#include "math.h"
#include "pxl/core.hh"
#include "pxl/hep.hh"


class Particle : public pxl::Particle {
  public:
    Particle();
    void setLabPosition(double, double, double, double);
    void setBjorkenX(double);
    void FullBoost(Particle*);
    void computeMomentumFromForce(pxl::Basic3Vector, double);
    void computeStep(double);

    double getRelativePtToCommonAxis(Particle*);
    double getDeltaR(Particle*);
    double getQ2(Particle*);

    double getBeta();
    double getGamma();

    double getX(size_t);
    double getLabTime();
    double getV(size_t);
    double getBjorkenX();
    
  private:
    std::vector<double> lab_x;
    double lab_t;
    double eigen_t;
    double bjorken_x;

};



#endif
