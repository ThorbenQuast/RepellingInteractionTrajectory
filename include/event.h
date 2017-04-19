#ifndef EVENT_H
#define EVENT_H

#include "math.h"
#include <cstdlib>
#include <ctime>

#include "particle.h"
#include "force.h"
#include "x_sec_weight.h"
#include "kinematics.h"
#include "constants.h"
#include "simResults.h"


class Event{
  public:
    Event();
    void setMaxIteration(double N_t);
    void setTMax(double tmax);
    void setMaxRadius(double maxR);
    void setEnergy(double E);
    void setAngleCut(double angleCut);
    void setCoincidenceTime(double coincidenceTime);
    void setRepulsion(double repulsion);

    void generate(simResult&);
    double computeStep(simResult&);
    double getWeight();

  private:
    double _tmax;
    int _N_t;
    double _maxRadius;
    double _E;
    double angleCut;
    double _coincidenceTime;
    bool repulsion;

    //book keeping for the simulation
    double weight;
    double _t;
    double _dt;
    int _t_counter;

    //generated objects
    Particle* ParticleA;
    Particle* ParticleB;
    double deltaT;

    xsecWeight* weights;
  
};


#endif