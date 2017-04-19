#include "force.h"
#include <iostream>

//pxl::Basic3Vector RelativisticForceLab(Particle* P, Particle* ParticleRef, bool repulsion) {
pxl::Basic3Vector Force(Particle* P, Particle* ParticleRef, bool repulsion) {
  //convert distances into natural unit!
  std::vector<double> dist;
  for (size_t i=0; i<3; i++) {
    dist.push_back((P->getX(i) - ParticleRef->getX(i))*meter_to_invGeV);
  }
  double magdist = sqrt(pow(dist[0],2)+pow(dist[1],2)+pow(dist[2],2));
  
  std::vector<double> beta;
  for (size_t i=0; i<3; i++) beta.push_back(ParticleRef->getV(i)/c);
  //for (size_t i=0; i<3; i++) beta.push_back(0.0);
  
  //formula for the strong force using a classical analogon
  double cF = repulsion ? +1./3. : -2./3.;  //colour factor and sign of the potential
  double alpha = 0.1190 * cF;  //natural units alpha_s=0.1190 is assumed  

  double commonPrefactor = alpha/pow(magdist-beta[0]*dist[0]-beta[1]*dist[1]-beta[2]*dist[2], 2);

  pxl::Basic3Vector F;
  F.setXYZ(commonPrefactor*(dist[0]/magdist-beta[0]), commonPrefactor*(dist[1]/magdist-beta[1]), commonPrefactor*(dist[2]/magdist-beta[2]));

  return F;
}


pxl::Basic3Vector ClassicForce(Particle* P, Particle* ParticleRef, bool repulsion) {
//pxl::Basic3Vector Force(Particle* P, Particle* ParticleRef, bool repulsion) {
  //convert distances into natural unit!
  std::vector<double> dist;
  for (size_t i=0; i<3; i++) dist.push_back((P->getX(i) - ParticleRef->getX(i))*meter_to_invGeV);
  double magdist = sqrt(pow(dist[0],2)+pow(dist[1],2)+pow(dist[2],2));
  

  //formula for the strong force using a classical analogon
  double cF = repulsion ? +1./3. : -2./3.;  //colour factor and sign of the potential
  double alpha = 0.1190 * cF;  //natural units alpha_s=0.1190 is assumed  

  double commonPrefactor = alpha/pow(magdist, 2);

  pxl::Basic3Vector F;
  F.setXYZ(commonPrefactor*(dist[0]/magdist), commonPrefactor*(dist[1]/magdist), commonPrefactor*(dist[2]/magdist));

  return F;
}


//deprecated
pxl::Basic3Vector _Force(Particle* P, Particle* ParticleRef, bool repulsion) {
//todo: debug the boost into rest frame

//compute forces there from the potential
  Particle* jet = new Particle();
  jet->setP4(P->getVector()+ParticleRef->getVector());

  //boost back into the lab frame (Minkowski force)
  Particle* boosted_P = new Particle();
  Particle* boosted_ParticleRef = new Particle();
  boosted_P->setP4(P->getVector());
  boosted_P->setLabPosition(P->getLabTime(), P->getX(0), P->getX(1), P->getX(2));
  boosted_ParticleRef->setP4(ParticleRef->getVector());
  boosted_P->setLabPosition(ParticleRef->getLabTime(), ParticleRef->getX(0), ParticleRef->getX(1), ParticleRef->getX(2));
  
  boosted_P->FullBoost(jet);
  boosted_ParticleRef->FullBoost(jet);

  //convert distances into natural unit!
  std::vector<double> dist;
  for (size_t i=0; i<3; i++) {
    dist.push_back((boosted_P->getX(i) - boosted_ParticleRef->getX(i))*meter_to_invGeV);
  }
  double magdist = sqrt(pow(dist[0],2)+pow(dist[1],2)+pow(dist[2],2));
  
  std::vector<double> beta;
  for (size_t i=0; i<3; i++) beta.push_back(boosted_ParticleRef->getV(i)/c);
  
  
  //formula for the strong force using a classical analogon
  double cF = repulsion ? +1./3. : -2./3.;  //colour factor and sign of the potential
  double alpha = 0.1190 * cF;  //natural units alpha_s=0.1190 is assumed  

  double commonPrefactor = alpha/pow(magdist-beta[0]*dist[0]-beta[1]*dist[1]-beta[2]*dist[2], 2);

  pxl::Basic3Vector F;
  F.setXYZ(commonPrefactor*(dist[0]/magdist-beta[0]), commonPrefactor*(dist[1]/magdist-beta[1]), commonPrefactor*(dist[2]/magdist-beta[2]));


  return F;
}