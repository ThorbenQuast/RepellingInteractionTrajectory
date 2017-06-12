#include "particle.h"
#include <iostream>

Particle::Particle(){
  lab_x.push_back(0); lab_x.push_back(0); lab_x.push_back(0); 
}

double Particle::getV(size_t i){
  double p_component = 0.;
  switch(i) {
    case 0:
      p_component = this->getPx();
      break;
    case 1:
      p_component = this->getPy();
      break;
    case 2:
      p_component = this->getPz();
      break;
  } 
  return p_component*c/this->getE();
}

void Particle::setLabPosition(double t, double x, double y, double z) {
  lab_t = t;
  lab_x[0] = x; lab_x[1] = y; lab_x[2] = z;
}

void::Particle::FullBoost(Particle* boost_particle) {
  this->boost(boost_particle->getBoostVector());

  std::cout<<boost_particle->getPx()<<std::endl;

  double beta = boost_particle->getBeta();
  double beta_x = boost_particle->getPx()/boost_particle->getE();
  double beta_y = boost_particle->getPy()/boost_particle->getE();
  double beta_z = boost_particle->getPz()/boost_particle->getE();
  double gamma = boost_particle->getGamma();

  std::cout<<"gamma: "<<gamma<<"   beta_x: "<<beta_x<<"   beta_y: "<<beta_y<<"   beta_z: "<<beta_z<<std::endl;

  double new_x0 = gamma*lab_t - beta_x*gamma*lab_x[0] - beta_y*lab_x[1] - beta_z*lab_x[2];
  double new_x1 = -beta_x*gamma*lab_t + (1+(gamma-1)*beta_x*beta_x/(beta*beta))*lab_x[0] + (gamma-1)*beta_x*beta_y/(beta*beta)*lab_x[1] + (gamma-1)*beta_x*beta_z/(beta*beta)*lab_x[2]; 
  double new_x2 = -beta_y*gamma*lab_t + (gamma-1)*beta_y*beta_x/(beta*beta)*lab_x[0] + (1.+(gamma-1)*beta_y*beta_y/(beta*beta))*lab_x[1] + (gamma-1)*beta_y*beta_z/(beta*beta)*lab_x[2]; 
  double new_x3 = -beta_z*gamma*lab_t + (gamma-1)*beta_z*beta_x/(beta*beta)*lab_x[0] + (gamma-1)*beta_z*beta_y/(beta*beta)*lab_x[1] + (1.+(gamma-1)*beta_z*beta_z/(beta*beta))*lab_x[2]; 
  
  std::cout<<"new_x1: "<<new_x1<<std::endl;
  std::cout<<"new_x2: "<<new_x2<<std::endl;
  std::cout<<"new_x3: "<<new_x3<<std::endl;
  lab_t = new_x0;
  lab_x[0] = new_x1;
  lab_x[1] = new_x2;
  lab_x[2] = new_x3;
}

void Particle::computeMomentumFromForce(pxl::Basic3Vector Force, double dt) {    //force is given in natural units
  double new_Px = this->getPx() + Force.getX()*dt*second_to_invGeV;
  double new_Py = this->getPy() + Force.getY()*dt*second_to_invGeV;
  double new_Pz = this->getPz() + Force.getZ()*dt*second_to_invGeV;
  double new_E = sqrt(pow(this->getMass(), 2) + pow(new_Px, 2) + pow(new_Py, 2) + pow(new_Pz, 2));

  this->setP4(new_Px, new_Py, new_Pz, new_E);

}

void Particle::setBjorkenX(double x) {
  bjorken_x = x;
}

double Particle::getBjorkenX() {
  return bjorken_x;
}

void Particle::computeStep(double dt) {
  lab_t += dt;
  std::cout<<"Lab t: "<<lab_t<<std::endl;
  eigen_t += dt/this->getGamma();
  double delta = 0;
  double p_component = 0;
  for (size_t i=0; i<3; i++) {
    switch(i) {
      case 0:
        p_component = this->getPx();
        break;
      case 1:
        p_component = this->getPy();
        break;
      case 2:
        p_component = this->getPz();
        break;
    }
    delta = c*p_component/this->getE()*dt;
    lab_x[i] = lab_x[i] + delta;
  }
}

double Particle::getX(size_t index) {
  return lab_x[index];
};

double Particle::getLabTime() {
  return lab_t;
}

double Particle::getBeta() {
  return this->getP()/this->getE();
}

double Particle::getGamma() {
  return this->getE()/this->getMass();
};

double Particle::getQ2(Particle* ref) {
  double Q2 = pow(this->getE()-ref->getE(), 2);
  Q2 -= pow(this->getPx()-ref->getPx(), 2);
  Q2 -= pow(this->getPy()-ref->getPy(), 2);
  Q2 -= pow(this->getPz()-ref->getPz(), 2);

  return -Q2;
}

double Particle::getRelativePtToCommonAxis(Particle* ref) {
  pxl::LorentzVector* main_axis = new pxl::LorentzVector(this->getVector()+ref->getVector());
  main_axis->normalize();
  double rel_pt = sqrt((pow(this->getP(),2)) - pow(main_axis->getX()*this->getPx() + main_axis->getY()*this->getPy() + main_axis->getZ()*this->getPz(),2));
  delete main_axis;
  return rel_pt;
}

double Particle::getDeltaR(Particle* ref) {
  double nominator = 0;
  for (size_t i=0; i<3; i++) {
    nominator+=ref->getV(i)*this->getV(i);
  }
  
  double len_v = 0;
  for (size_t i=0; i<3; i++) len_v+=this->getV(i)*this->getV(i);
  len_v = sqrt(len_v);

  double len_ref = 0;
  for (size_t i=0; i<3; i++) len_ref+=ref->getV(i)*ref->getV(i);
  len_ref = sqrt(len_ref);
  
  return acos(nominator/(len_v*len_ref));
}

