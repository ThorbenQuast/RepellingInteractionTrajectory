#include "x_sec_weight.h"
#include <iostream>


//chosen parameterization of the parton density function
double a(double x) {
  double logx10 = log(x)/log(10.0);
  return 0.1 + p0_pdf_a * exp(-p1_pdf_a * ((p2_pdf_a-logx10)+ exp(-(p2_pdf_a - logx10)))); 
};

double kappa(double x) {
  double logx10 = log(x)/log(10.0);
  return p0_pdf_kappa + p1_pdf_kappa * pow(((1-logx10)*p2_pdf_kappa), p3_pdf_kappa);  
}

double xpdf(double x, double t) {
  return a(x)*pow(log(t/l_pdf),kappa(x));
};
  

xsecWeight::xsecWeight() {
  m_source = 0;
  m_target = 0;
  E_source = 0;
  s = 0;
  p = 0;
  _beta = 0;
  _gamma = 0;
  boosted_Ep = 0;
  boosted_Ee = 0;
};

void xsecWeight::setSourceMass(double mass) {
  m_source = mass;
  _updateParams();
};

void xsecWeight::setTargetMass(double mass) {
  m_target = mass;
  _updateParams();
}; 

void xsecWeight::setSourceEnergy(double energy) {
  E_source = energy;
  _updateParams();
};

void xsecWeight::_updateParams(){
  try { 
    s = 2*m_target*E_source+pow(m_target,2);
    p = pow((pow(E_source,2)-pow(m_source,2)),0.5);
    _beta = (p)/(m_target+E_source);
    _gamma = 1.0/sqrt(1-pow(_beta,2));
    boosted_Ep = _gamma*(m_target);
    boosted_Ee = _gamma*(E_source-_beta*p);
  }
  catch(int e) {
    return;
  }
};

double xsecWeight::_x(double theta_q, double t ) {
  double eta = eta_q(theta_q);
  return 1.0/(2*boosted_Ep*s) * (boosted_Ep*t + sqrt(boosted_Ep*t*(boosted_Ep*t+4.0*s*boosted_Ee*exp(2.0*eta))));
};

double xsecWeight::_d2sigma_dthetaq_dt(double theta_q, double t) {
  if (t < 3.5) return 0.0;
  double x = _x(theta_q, t);
  if (!(x<1.0 and x>0)) {
    std::cout<<"x is greater than 0"<<std::endl;
    return 0.0;
  } 

  double eta = eta_q(theta_q);
  double main = (pow((1./137.),2))/(8.0*pi) * 1.0/(x*(pow(t,2)))*xpdf(x, t); 
  double jacobian = 2.0*exp(2.0*eta)*t*boosted_Ee/sqrt(boosted_Ep*t*(boosted_Ep*t+4.0*s*boosted_Ee*exp(2.0*eta))) * 1.0/sin(theta_q);
  return main*jacobian*0.389*((10.0e-31));
};

