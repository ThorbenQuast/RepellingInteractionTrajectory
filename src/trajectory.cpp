#include "trajectory.h"

std::pair<double, double> runTrajectorySimulation(bool repulsion, double Energy, double radialDistance, double coincidenceTime, double angleCut, double electronBeamAngle, double T_max, simResult &results) {
  Event* event = new Event();

  //set the experimental conditions
  event->setEnergy(Energy);
  event->setMaxRadius(radialDistance);
  event->setCoincidenceTime(coincidenceTime);
  event->setAngleCut(angleCut);
  event->setElectronBeamAngle(electronBeamAngle);

  //generate it
  event->generate(results);

  //settings for the simulation of the trajectory
  int maxIteration=1000;
  event->setMaxIteration(maxIteration);  //maximum nummber of iterations before steps in dt are increased by x10
  results.set_maxIteration(maxIteration);
  event->setTMax(T_max); //#this corresponds to a flight distance of ~3m, should be long enough
  event->setRepulsion(repulsion);

  //run the simulation
  double status;
  do {
    status = event->computeStep(results);
  } while(status==-1);

  delete event;

  return std::make_pair(status, event->getWeight());
};
