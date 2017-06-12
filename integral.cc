//compile and run:
//g++ integral.cc -o integral.out && ./integral.out

#include <iostream>
#include <cmath>

double c = 299792458.0; //speed of light in m/s
double m_e = 0.511*10e-3; //electron mass in keV

double pi = 3.1415926535897;

double deg_to_rad = (pi/180);
double invGev_to_m2 = 0.389379*pow(10, -3) * pow(10, -28);

/**********/
//parameterisation of the structure function
double p0_pdf_a = 1.42334682e+03;
double p1_pdf_a = 8.13911499e+00;
double p2_pdf_a = -6.44188039e-01;

double p0_pdf_kappa = -5.27385445e+02;
double p1_pdf_kappa = 5.24698650e+02;
double p2_pdf_kappa = 1.98072582e+00;
double p3_pdf_kappa = 3.66260610e-03;
double l_pdf = 0.35;


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
/**********/


double eta_q(double theta) {
  return -log(tan(theta/2));
};  


//cross section formula
double d2sigma_dthetaq_dt(double bjorken_x, double t) {  
	if (bjorken_x <= 0 || bjorken_x > 1) return 0.0;
  if (t < 3.5) return 0.0;

  return xpdf(bjorken_x, t)/(bjorken_x * pow(t,2));
};


//checked with Christian's "get_x" function --> same result
double get_x(double eta_q, double Q2, double s, double Ep, double Ee) {
	return Q2/(2.*s) + sqrt(pow(Q2/(2.*s), 2) + Q2/s * Ee/Ep * exp(2*eta_q));
}


/**********/
//differentials to be multiplied to the integrand
double get_dx_detaq(double eta_q, double Q2, double s, double Ep, double Ee) {
	return Q2/s * Ee/Ep * exp(2*eta_q) / sqrt(pow(Q2/(2*s), 2) + Q2/s * Ee/Ep * exp(2*eta_q));
} //cross-checked with Wolfram alpha

double get_detaq_dx(double theta_q) {
	return -1./tan(theta_q/2.) / pow(cos(theta_q/2.), 2) / 2.;			//note: equivalent to -/sin(theta)
}
/**********/



int main(int argc, char* argv[]){
	//machine parameters
	double Ep = 1.0;		//proton's mass is 1GeV
	double Ee = 60.;		//electron energy is set to 60GeV
	double s = 2.*Ee*Ep + pow(Ep, 2);	//center-of-mass energy squared

	int N_granularity = 20000;		//number of steps, adjust here.

	double Q2, Q2_start; Q2 = Q2_start = 3.5;		//minimum momentum transfer in DIS					
	double theta_q, theta_q_start; theta_q = theta_q_start = 0.01*deg_to_rad;		//no scattering into the beam-pipe

	double dQ2 = (pow(Ee, 2) - Q2)/N_granularity;
	double dtheta_q = (180.*deg_to_rad-2*theta_q_start) / N_granularity;
  

	//integrate as function of theta and Q2
	double integral = 0.;

	//test cout:
	Q2 = 10.;
	theta_q = 2.6;
	double _eta_q_test = eta_q(theta_q);
	double x_test = get_x(_eta_q_test, Q2, s, Ep, Ee);
	std::cout<<"Comparison with Christian: Q2 = "<< Q2<<" theta_q = "<<theta_q<<" --> "<<d2sigma_dthetaq_dt(x_test, Q2) * fabs(get_dx_detaq(_eta_q_test, Q2, s, Ep, Ee)) * fabs(get_detaq_dx(theta_q)) * pow(1./137, 2) * (1./(4.*pi)) * invGev_to_m2*pow(10, 28)<<std::endl;

	
	theta_q = theta_q_start;
	//integral is approximated as a sum of all the rectangles
	for (int i=0; i<N_granularity; i++) {
		theta_q = theta_q + dtheta_q;
		Q2 = Q2_start;
		for (int j=0; j<N_granularity; j++) {
			Q2 = Q2 + dQ2;
			
			double _eta_q = eta_q(theta_q);
			double x = get_x(_eta_q, Q2, s, Ep, Ee);
			integral = integral + d2sigma_dthetaq_dt(x, Q2) * fabs(get_dx_detaq(_eta_q, Q2, s, Ep, Ee)) * fabs(get_detaq_dx(theta_q)) * dtheta_q * dQ2;
		}
	}


	integral = pow(1./137, 2) * (1./(4.*pi)) * integral * invGev_to_m2;
	
	std::cout<<"Integrated cross section: "<<integral<<" m2 = "<<integral*pow(10, 28+12)<< "pb" <<std::endl;

	
	return 0;
}
