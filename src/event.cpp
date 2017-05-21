#include "event.h"
#include <iostream>

Event::Event() {
    _tmax = 0;
    _N_t = 0;
    _maxRadius = 0;
    _E = 0;
    _s = 0.;
    _electronBeamAngle = 0.;    //mixing angle between y and z momentum components
    angleCut = pi;    
    _coincidenceTime = 0.;  
    repulsion = false;  

    //book keeping for the simulation
    weight = 0.;
    _t = 0.;
    _dt = 1e-34;   
    _t_counter = 0;

    //generated objects
    ParticleA = new Particle();
    ParticleB = new Particle();
    deltaT = 0.;    //time in between two collisions

};


void Event::setMaxIteration(double N_t) {
  _N_t = N_t;
}

void Event::setTMax(double tmax) {
  _tmax = tmax;
}

void Event::setMaxRadius(double maxR) {
  _maxRadius = maxR;
}

void Event::setEnergy(double E) {
  _E = E;
  _s = 2 * 0.931 * E;   //proton mass is assumed
}

void Event::setElectronBeamAngle(double beamAngle) {
  _electronBeamAngle = beamAngle*pi/180.;
}

void Event::setAngleCut(double angleCutDegree) {  //argument is in degrees
  angleCut = angleCutDegree*pi/180.;
}

void Event::setCoincidenceTime(double coincidenceTime) {
  _coincidenceTime = coincidenceTime;
}

void Event::setRepulsion(double rep) {
  repulsion = rep;
};

double Event::getWeight() {
  return weight;
}

void Event::generate(simResult &results) {
  
  //define the quark masses to be other from the up or bottom
  double randMassA = rand()*1./(RAND_MAX);
  double massAGeV = randMassA<0.5 ? 0.0023 : 0.0048;

  double randMassB = rand()*1./(RAND_MAX);
  double massBGeV = randMassB<0.5 ? 0.0023 : 0.0048;


  //generate the relative displacements, ParticleA is always defined to be in the origin
  ParticleA->setLabPosition(0.0, 0.0, 0.0, 0.0);

  double r_norm, r_try;
  do {
   r_norm = rand()*1./(RAND_MAX);
   r_try = rand()*1./(RAND_MAX);
  } while(r_norm < r_try);

  double phi_displacement = rand()*2.*pi/(RAND_MAX);

  double theta_displacement = acos(rand()*2./(RAND_MAX)-1.);   //from : http://mathworld.wolfram.com/SpherePointPicking.html


  double deltaX12 = _maxRadius*r_norm*cos(phi_displacement)*sin(theta_displacement);
  double deltaY12 = _maxRadius*r_norm*sin(phi_displacement)*sin(theta_displacement);
  double deltaZ12 = _maxRadius*r_norm*cos(theta_displacement);

  ParticleB->setLabPosition(0.0, deltaX12, deltaY12, deltaZ12);

  //set the coincidence time:
  deltaT = rand()*1./(RAND_MAX)*_coincidenceTime;


  double p_electron = sqrt(pow(_E, 2) - pow(511000.0e-9, 2));
  
  /**********************/
  //generate the kinematics!
  Particle* IncidentElectron1 = new Particle();
  IncidentElectron1->setP4(0., 0., p_electron, _E);
  Particle* OriginalQuarkA = new Particle();
  double px_A_orig = 0.2*(rand()*1./(RAND_MAX)-0.5)/sqrt(3.);
  double py_A_orig = 0.2*(rand()*1./(RAND_MAX)-0.5)/sqrt(3.);
  double pz_A_orig = 0.2*(rand()*1./(RAND_MAX)-0.5)/sqrt(3.);
  double E_A_orig = sqrt(pow(massAGeV,2)+pow(px_A_orig,2)+pow(py_A_orig,2)+pow(pz_A_orig,2));
  OriginalQuarkA->setP4(px_A_orig, py_A_orig, pz_A_orig, E_A_orig);

  Particle* outgoingElectron1 = new Particle();
  generateKinematic(IncidentElectron1, OriginalQuarkA, outgoingElectron1, ParticleA);
  //outgoingElectron1->setP4(-1.01304, 1.35215, 43.8244, 43.8569);
  //ParticleA->setP4(1.01207, -1.29981, 16.1292, 16.2131);
  //ParticleA->setP4(0., 0., 3.5, sqrt(3.5*3.5+pow(massAGeV,2)));

  double QA = outgoingElectron1->getQ2(IncidentElectron1);
  double thetaA = ParticleA->getTheta();


  Particle* IncidentElectron2 = new Particle();
  IncidentElectron2->setP4(0., sin(_electronBeamAngle)*p_electron, cos(_electronBeamAngle)*p_electron, _E);
  
  //dice the kinematics of the second particle according to its electron
  Particle* OriginalQuarkB = new Particle();
  double px_B_orig = 0.2*(rand()*1./(RAND_MAX)-0.5)/sqrt(3.);
  double py_B_orig = 0.2*(rand()*1./(RAND_MAX)-0.5)/sqrt(3.);
  double pz_B_orig = 0.2*(rand()*1./(RAND_MAX)-0.5)/sqrt(3.);
  double E_B_orig = sqrt(pow(massBGeV,2)+pow(px_B_orig,2)+pow(py_B_orig,2)+pow(pz_B_orig,2));
  OriginalQuarkB->setP4(px_B_orig, py_B_orig, pz_B_orig, E_B_orig);

  Particle* outgoingElectron2 = new Particle();
  generateKinematic(IncidentElectron2, OriginalQuarkB, outgoingElectron2, ParticleB);
  //outgoingElectron2->setP4(-0.729088, 1.19555, 31.0232, 31.0548 );
  //ParticleB->setP4(0.780043, -1.10741, 29.0273, 29.0589);
  //ParticleB->setP4(0., 0., 3.5, sqrt(3.5*3.5+pow(massBGeV,2)));
  double QB = outgoingElectron2->getQ2(IncidentElectron2);
  double thetaB = ParticleB->getTheta();

  /**********************/


  weight = d2sigma_dthetaq_dt(OriginalQuarkA->getBjorkenX(), QA, _s) * d2sigma_dthetaq_dt(OriginalQuarkB->getBjorkenX(), QB, _s);

  //cut on outgoing quark flight angle w.r.t. each other
  if (ParticleA->getDeltaR(ParticleB) > angleCut) weight = 0.;


  if (weight>0) {
    std::cout<<"DeltaR = "<<outgoingElectron1->getDeltaR(outgoingElectron2)<<std::endl;
    std::cout<<"Relative pt at the beginning: "<<ParticleA->getRelativePtToCommonAxis(ParticleB)<<std::endl;
    std::cout<<"INIT"<<std::endl;
    std::cout<<"Particle 1"<<std::endl;
    std::cout<<"electron E: "<<outgoingElectron1->getE()<<"    quark E: "<<ParticleA->getE()<<std::endl;
    std::cout<<"electron px: "<<outgoingElectron1->getPx()<<"    quark px: "<<ParticleA->getPx()<<"    quark x: "<<ParticleA->getX(0)<<std::endl;
    std::cout<<"electron py: "<<outgoingElectron1->getPy()<<"    quark py: "<<ParticleA->getPy()<<"    quark y: "<<ParticleA->getX(1)<<std::endl;
    std::cout<<"electron pz: "<<outgoingElectron1->getPz()<<"    quark pz: "<<ParticleA->getPz()<<"    quark z: "<<ParticleA->getX(2)<<std::endl;
    std::cout<<"QA: "<<QA<<"   thetaA:  "<<thetaA<<std::endl;

    std::cout<<"Particle 2"<<std::endl;
    std::cout<<"electron E: "<<outgoingElectron2->getE()<<"    quark E: "<<ParticleB->getE()<<std::endl;
    std::cout<<"electron px: "<<outgoingElectron2->getPx()<<"    quark px: "<<ParticleB->getPx()<<"    quark x: "<<ParticleB->getX(0)<<std::endl;
    std::cout<<"electron py: "<<outgoingElectron2->getPy()<<"    quark py: "<<ParticleB->getPy()<<"    quark y: "<<ParticleB->getX(1)<<std::endl;
    std::cout<<"electron pz: "<<outgoingElectron2->getPz()<<"    quark pz: "<<ParticleB->getPz()<<"    quark z: "<<ParticleB->getX(2)<<std::endl;
    std::cout<<"QB: "<<QB<<"   thetaB:  "<<thetaB<<std::endl;
    std::cout<<std::endl;
  
    results.addEntry("mass_q1", massAGeV);
    results.addEntry("mass_q2", massBGeV);
    results.addEntry("IP1_x", 0.);
    results.addEntry("IP1_y", 0.);
    results.addEntry("IP1_z", 0.);
    results.addEntry("IP2_x", deltaX12);
    results.addEntry("IP2_y", deltaY12);
    results.addEntry("IP2_z", deltaZ12);
    results.addEntry("px_e1_in", IncidentElectron1->getPx());
    results.addEntry("py_e1_in", IncidentElectron1->getPy());
    results.addEntry("pz_e1_in", IncidentElectron1->getPz());
    results.addEntry("E_e1_in", IncidentElectron1->getE());
    results.addEntry("px_e2_in", IncidentElectron2->getPx());
    results.addEntry("py_e2_in", IncidentElectron2->getPy());
    results.addEntry("pz_e2_in", IncidentElectron2->getPz());
    results.addEntry("E_e2_in", IncidentElectron2->getE());
    results.addEntry("px_q1_in", px_A_orig);
    results.addEntry("py_q1_in", py_A_orig);
    results.addEntry("pz_q1_in", pz_A_orig);
    results.addEntry("E_q1_in", E_A_orig);
    results.addEntry("x_q1", OriginalQuarkA->getBjorkenX());
    results.addEntry("px_e1_out", outgoingElectron1->getPx());
    results.addEntry("py_e1_out", outgoingElectron1->getPy());
    results.addEntry("pz_e1_out", outgoingElectron1->getPz());
    results.addEntry("E_e1_out", outgoingElectron1->getE());
    results.addEntry("px_q1_out", ParticleA->getPx());
    results.addEntry("py_q1_out", ParticleA->getPy());
    results.addEntry("pz_q1_out", ParticleA->getPz());
    results.addEntry("E_q1_out", ParticleA->getE());
    results.addEntry("Q2_1", QA);
    results.addEntry("px_q2_in", px_B_orig);
    results.addEntry("py_q2_in", py_B_orig);
    results.addEntry("pz_q2_in", pz_B_orig);
    results.addEntry("E_q2_in", E_B_orig);
    results.addEntry("x_q2", OriginalQuarkB->getBjorkenX());
    results.addEntry("px_e2_out", outgoingElectron2->getPx());
    results.addEntry("py_e2_out", outgoingElectron2->getPy());
    results.addEntry("pz_e2_out", outgoingElectron2->getPz());
    results.addEntry("E_e2_out", outgoingElectron2->getE());
    results.addEntry("px_q2_out", ParticleB->getPx());
    results.addEntry("py_q2_out", ParticleB->getPy());
    results.addEntry("pz_q2_out", ParticleB->getPz());
    results.addEntry("E_q2_out", ParticleB->getE());
    results.addEntry("Q2_2", QB);
    results.addEntry("DeltaR_initial", ParticleA->getDeltaR(ParticleB));
    results.addEntry("relPt_initial", ParticleA->getRelativePtToCommonAxis(ParticleB));
    results.addEntry("event_weight", weight);
  }

  delete IncidentElectron1;
  delete IncidentElectron2;
  delete outgoingElectron1;
  delete outgoingElectron2;
}


double Event::computeStep(simResult &results) {
  //if (_t > _tmax) return ParticleA->getDeltaR(ParticleB); //return the final value
  if (_t > _tmax) {
    std::cout<<"Final kinematics: "<<std::endl;
    std::cout<<"Particle 1"<<std::endl;
    std::cout<<"    quark px: "<<ParticleA->getPx()<<std::endl;
    std::cout<<"    quark py: "<<ParticleA->getPy()<<std::endl;
    std::cout<<"    quark pz: "<<ParticleA->getPz()<<std::endl;
    std::cout<<"Particle 2"<<std::endl;
    std::cout<<"    quark px: "<<ParticleB->getPx()<<std::endl;
    std::cout<<"    quark py: "<<ParticleB->getPy()<<std::endl;
    std::cout<<"    quark pz: "<<ParticleB->getPz()<<std::endl;

    std::cout<<"rel pt 1: "<<ParticleA->getRelativePtToCommonAxis(ParticleB)<<std::endl;
    std::cout<<"rel pt 2: "<<ParticleB->getRelativePtToCommonAxis(ParticleA)<<std::endl;
    
  
    results.addEntry("px_q1_final", ParticleA->getPx());
    results.addEntry("py_q1_final", ParticleA->getPy());
    results.addEntry("pz_q1_final", ParticleA->getPz());
    results.addEntry("E_q1_final", ParticleA->getE());
    results.addEntry("px_q2_final", ParticleB->getPx());
    results.addEntry("py_q2_final", ParticleB->getPy());
    results.addEntry("pz_q2_final", ParticleB->getPz());
    results.addEntry("E_q2_final", ParticleB->getE());
    results.addEntry("relPt_final", ParticleA->getRelativePtToCommonAxis(ParticleB));
    
    return ParticleA->getRelativePtToCommonAxis(ParticleB); //return the final value
  }
  //abort condition: in case the event is to be discarded <-->has weight 0
  // - either the kinematic situation is not allowed (see x-sec calculation)
  // - the angular separation in the flight direction is too high from the beginning (indicated parameter)
  if (weight==0) return -999;   //abort condition

  //only propagate first particle freely if the second one is not generated yet
  if (_t < _coincidenceTime)
    ParticleA->computeStep(_dt);
  else {
    pxl::Basic3Vector ForceA = Force(ParticleA, ParticleB, repulsion);
    pxl::Basic3Vector ForceB = Force(ParticleB, ParticleA, repulsion);

    ParticleA->computeMomentumFromForce(ForceA, _dt);
    ParticleB->computeMomentumFromForce(ForceB, _dt);

    ParticleA->computeStep(_dt);
    ParticleB->computeStep(_dt);
  }

  _t += _dt;
  _t_counter += 1;

  //the stepsize in (laboratory) time varies as more precision is needed at close distances of the quarks to each other (huge forces->accelerations)
  //therefore, dt is adjusted by x10 every _N_t simulation steps
  if (!(_t_counter % _N_t)) { //then increase the stepsize dt
    _dt*=10.;
  }

  return -1;  //status for computation
};

