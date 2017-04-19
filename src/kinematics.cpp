#include "kinematics.h"

//computes the angles of the scattering process respecting conservation of 4-momentum

//a=sum py, b=sum px, X=pt electron, Y = pt quark
//phi=angle electron, theta=angle quark
std::pair<double, double> compute_phi_theta(double a, double b, double X, double Y) {
  double theta;double phi;

  double left_side = (pow(a, 2)+pow(b, 2)-pow(X, 2)+pow(Y, 2))/(2*Y*sqrt(pow(a, 2)+pow(b, 2)));
  if (a>0)
    theta = asin(left_side)-atan(b/a);
  else if (b>0)
    theta = acos(left_side)+atan(a/b);
  else {
    theta = asin(left_side)-(atan(b/a)-pi);
  }


  double sin_phi = (a-Y*sin(theta))/X;
  double cos_phi = (b-Y*cos(theta))/X;


  if (sin_phi>0 && cos_phi>0)
    phi = asin(sin_phi);
  else if (sin_phi>0)
    phi = acos(cos_phi);
  else if (cos_phi>0)
    phi = asin(sin_phi);
  else 
    phi = asin(sin_phi)+pi/2;


  return std::make_pair(phi, theta);
}


void generateKinematic(Particle* IncidentElectron, Particle* OriginalQuark, Particle* outgoingElectron, Particle* Quark) {
  double sum_E = IncidentElectron->getE()+OriginalQuark->getE();
  double sum_px = IncidentElectron->getPx()+OriginalQuark->getPx();
  double sum_py = IncidentElectron->getPy()+OriginalQuark->getPy();
  double sum_pz = IncidentElectron->getPz()+OriginalQuark->getPz();

  double fraction_Energy;
  double fraction_Pz;

  double electron_E;
  double quark_E;
  double electron_pz;
  double quark_pz;

  double electron_pt;
  double quark_pt;

  bool invalidKinematics = true;
  do {
    fraction_Energy = rand()*1./(RAND_MAX);
    fraction_Pz = rand()*1./(RAND_MAX);

    electron_E = fraction_Energy*sum_E;
    quark_E = (1-fraction_Energy)*sum_E;
    electron_pz = fraction_Pz*sum_pz;
    quark_pz = (1-fraction_Pz)*sum_pz;

    electron_pt = sqrt(pow(electron_E,2)-pow(outgoingElectron->getMass(),2)-pow(electron_pz,2));
    quark_pt = sqrt(pow(quark_E,2)-pow(Quark->getMass(),2)-pow(quark_pz,2));
  
    invalidKinematics = (sqrt(pow(electron_E,2)-pow(outgoingElectron->getMass(),2))<electron_pz);
    invalidKinematics = invalidKinematics || (sqrt(pow(quark_E,2)-pow(Quark->getMass(),2))<quark_pz);
    invalidKinematics = invalidKinematics || (pow(sum_px,2)+pow(sum_py,2) > pow(electron_pt+quark_pt,2));
    invalidKinematics = invalidKinematics || (pow(sum_px,2)+pow(sum_py,2) < pow(electron_pt-quark_pt,2));
  } while(invalidKinematics);

  //determine the transverse kinematics
  
  std::pair<double, double> electronAngle_quarkAngle = compute_phi_theta(sum_py, sum_px, electron_pt, quark_pt);
  double electron_px = electron_pt*cos(electronAngle_quarkAngle.first);
  double electron_py = electron_pt*sin(electronAngle_quarkAngle.first);
  double quark_px = quark_pt*cos(electronAngle_quarkAngle.second);
  double quark_py = quark_pt*sin(electronAngle_quarkAngle.second);

  //std::cout<<"sum_px: "<<sum_px<<"   electron_px: "<<electron_px<<"   quark px: "<<quark_px<<std::endl;
  //std::cout<<"sum_py: "<<sum_py<<"   electron_py: "<<electron_py<<"   quark py: "<<quark_py<<std::endl;
  //std::cout<<"sum_pz: "<<sum_pz<<"   electron_pz: "<<electron_pz<<"   quark pz: "<<quark_pz<<std::endl;
  //std::cout<<"sum_E: "<<sum_E<<"   electron_E: "<<electron_E<<"   quark E: "<<quark_E<<std::endl;

  outgoingElectron->setP4(electron_px, electron_py, electron_pz, electron_E);
  Quark->setP4(quark_px, quark_py, quark_pz, quark_E);

};

double eta_q(double theta_q) {
    return -log(tan(theta_q/2.0));
};

