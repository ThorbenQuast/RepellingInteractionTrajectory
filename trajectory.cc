//main
#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include "trajectory.h"

#include "StrongInteractionTrajectoryConfig.h"

#include "simResults.h"

int main(int argc, char* argv[]){
  std::cout<<"StrongInteractionTrajectory v."<<StrongInteractionTrajectory_VERSION_MAJOR<<"."<<StrongInteractionTrajectory_VERSION_MINOR<<std::endl;
  if (argc < 2){
    std::cout<<"Correct usage: OutputFile(<path>) N_sim(<integer>)"<<std::endl;
    return 1;
  }
  bool repulsion = true;
  int N_sim = atoi(argv[2]);

  std::cout<<"Running with repulsion: "<<repulsion<<"   creating "<<argv[1]<<"  and simulating "<<N_sim<<" trajectories"<<std::endl;

  //intialize the random number generator
  srand (time(NULL));

  //configuration
  double Energy = 60.;      //in GeV, energy of incoming electron beam
  double radialDistance = 1e-15;   //maximum distance of the two generated quarks in the plane perpendicular to the beam axis
  //double coincidenceTime = 10.0e-23; //equivalent to the time bin = coincidence time of both scattering processes, maximum time interval in which the two quarks are creared
  double coincidenceTime = 0.; //equivalent to the time bin = coincidence time of both scattering processes, maximum time interval in which the two quarks are creared
  double angleCut = 10.;    //initial direction angle cut (if the angle between the quarks is higher after generation, trajectories are not computed)
  double electronBeamAngle = 1.;  //angle between the incident electron beams
  double T_max = 1e-8;   //maximum time to simulate the trajectory


  simResult results;
  results.set_energy(Energy);
  results.set_radialDistance(radialDistance);
  results.set_coincidenceTime(coincidenceTime);
  results.set_electronBeamAngle(electronBeamAngle);
  results.set_angleCut(angleCut);
  results.set_T_max(T_max);
  results.set_repulsion(repulsion);

  std::pair<double, double> _event_result;
  for (size_t i=0; i<N_sim; ) {
    _event_result = runTrajectorySimulation(repulsion, Energy, radialDistance, coincidenceTime, angleCut, electronBeamAngle, T_max, results);
    if (_event_result.first >= 0) {
      i++;
      std::cout<<std::endl<<"  ****  "<<std::endl<<std::endl;;
      std::cout<<"Iteration "<<i<<":    relPt="<<_event_result.first<<"   weight="<<_event_result.second<<std::endl;
      std::cout<<std::endl<<"  ****  "<<std::endl;
    }
  }

  std::ofstream outfile;
  outfile.open(argv[1], std::ofstream::out | std::ofstream::trunc);
  results.printHeader(outfile);
  results.printEntries(outfile);
  outfile.close();
  return 0;
}
