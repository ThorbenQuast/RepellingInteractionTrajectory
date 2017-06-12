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
  return a(x)*pow(2*log(t/l_pdf),kappa(x));
};

double eta_q(double theta) {
  return -log(tan(theta/2));
};  


double d2sigma_dthetaq_dt(double bjorken_x, double t) {  
	if (bjorken_x <= 0 || bjorken_x > 1) return 0.0;
  if (t < 3.5) return 0.0;
  return xpdf(bjorken_x, t)/(bjorken_x * pow(t,2));
};
