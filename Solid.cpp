// Solid.cpp
//
// 
//This code computes the solid angle subtended by a 
// rapidly rotating star.
//
//
// Based on code written by Coire Cadeau.
//
// (C) Coire Cadeau, 2007

// Source (C) Coire Cadeau 2007, all rights reserved.
//
// Permission is granted for private use only, and not
// distribution, either verbatim or of derivative works,
// in whole or in part.
//
// The code is not thoroughly tested or guaranteed for
// any particular use.

// Full spectral parameters and 2 energy bands included.

#include <iostream>
#include <fstream>
#include <cmath>
#include <exception>
#include <vector>
#include <string>
#include "OblDeflectionTOA.h"
#include "Chi.h"
#include "PolyOblModelNHQS.h"
#include "PolyOblModelCFLQS.h"
#include "SphericalOblModel.h"
#include "OblModelBase.h"
#include "Units.h"
#include "Exception.h"
#include "Struct.h"




int main(int argc, char** argv) try {

  int nmu(400);
 
  double incl, theta(90.0), mass, rspot(0.0), omega, req, bbrat, ts;
  double distance, Temperature, E0(1.0);
  double north, south;
  double omega_cgs;
  unsigned int modelchoice(1), quadmodel(0), fdmodel(0);
  int run;
  unsigned int numbins(0);
  bool incl_is_set(false), theta_is_set(false), bbrat_is_set(false), numbins_is_set(false),
    mass_is_set(false), rspot_is_set(false), omega_is_set(false),
    model_is_set(true), infile_is_set(false), ignore_time_delays(false);
  char out_file[256] = "oblflux.txt";

  int power=0;

  char data_file[80];
  char line[80];
  char vert[180];
  char min_file[180];
  //char min_store[180];
  char out_dir[80];
  // New variables added by SMM
 
  double aniso(0.586), Gamma(2.0);
  double shift_t(0.0);    /* Time off-set from data */
  double ftol = 1e-3;
  double b_mid;
  std::ifstream data;
  std::ofstream out;

  // Create LightCurve data structure
  class LightCurve curve;
  double mu, cosg;
 
  curve.para_read_in = false;


  // Read in the command line options
  for(int i(1); i < argc; i++) {
    if(argv[i][0]=='-') {
      switch(argv[i][1]) {

      case 'o':
	// Name of output file
	sscanf(argv[i+1], "%s", out_file);
	break;

      case 's':
	// This option controls the use of an external file
	// to define the shape of the spot.
	infile_is_set = true;
	break;

      case 'q':
	// Oblateness model (default is 1)
	sscanf(argv[i+1], "%u", &modelchoice);
	model_is_set = true;
	break;

      case 'Q':
	// Quadrupole model (default is 0)
	sscanf(argv[i+1], "%u", &quadmodel);
	break;

      case 'w':
	// Frame dragging (default is 0)
	sscanf(argv[i+1],"%u",&fdmodel);
	break;

      case 'n':
	//number of bins
	sscanf(argv[i+1], "%u", &numbins);
	numbins_is_set = true;
	break;

      case 'p':
	//power that eta is raised to
	sscanf(argv[i+1], "%d", &power);
	break;

      case 'i':
	// Inclination angle of the observer (degrees)
	sscanf(argv[i+1], "%lf", &incl);
	incl_is_set = true;
	break;

      case 'e':
	// Emission angle of the spot (degrees)
	sscanf(argv[i+1], "%lf", &theta);
	theta_is_set = true;
	break;

      case 'm':
	// Mass of the star (solar mass units)
	sscanf(argv[i+1], "%lf", &mass);
	mass_is_set = true;
	break;

      case 'r':
	// Radius of the star at the equator(km)
	sscanf(argv[i+1], "%lf", &req);
	rspot_is_set = true;
	break;

      case 'f': 
	// Spin frequency (Hz)
	sscanf(argv[i+1], "%lf", &omega);
	omega_cgs = omega;
	omega_is_set = true;
	break;

      case 't': //toggle ignore_time_delays (only affects output)
	ignore_time_delays = true;
	break;

      case 'T':  // Temperature of the spot, in the star's frame, in keV
	sscanf(argv[i+1], "%lf", &Temperature);
	break;

      case 'D':  // Distance to NS
	sscanf(argv[i+1], "%lf", &distance);
	break;

      case 'h':
      
	std::cout << "solid help: " << std::endl
		  << "-o filename: (optional, default is oblflux.txt)" << std::endl
		  << "-s filename: (optional)" <<std::endl
		  << "-n numbins: number of datapoints" << std::endl
		  << "-i inclination of observer in degrees, between 0 and 180." << std::endl
		  << "-e location of emission region in degrees, between 0 and 180." << std::endl
		  << "-m mass of star in Msun." << std::endl
		  << "-r radius of star (at the spot) in km." << std::endl
		  << "-f rotation frequency of star in Hz." << std::endl
		  << "-t ignore time delays in output (see source)." << std::endl
		  << "-p azimuthal location (phi_0) of spot [0.0]." << std::endl
		  << "-q #, where #=:" << std::endl
		  << "      1 for Neutron/Hybrid quark star poly model" << std::endl
		  << "      2 for CFL quark star poly model" << std::endl
		  << "      3 for spherical model" << std::endl
		  << "-Q [0] Quadrupole model 0 or 1" << std::endl
		  << std::endl;
	return 0;
      } // end switch	
    } // end if
  } // end for

  // Check that necessary numbers are set.
  if( !( incl_is_set 
	 && numbins_is_set
	 && mass_is_set
	 && rspot_is_set
	 && omega_is_set
	 ) ) {
    std::cout << "Not all required parameters were specified. Exiting."
	      << std::endl;
    return -1;
  }

  // Print out information about the model
  std::cout << std::endl << "Pulse: Stellar Model Parameters" << std::endl << std::endl;
  std::cout << "Mass = " << mass << " Msun" << std::endl;
  std::cout << "Req = " << req << " km" <<  std::endl;
  std::cout << "2GM/Rc^2 = " << 2.0*Units::G*mass*Units::MSUN/(req*1e+5*Units::C*Units::C)  <<  std::endl;
  std::cout << "spin = " << omega << " Hz" << std::endl;
  std::cout << "v/c = " << req * 1e5 * omega * 2.0*Units::PI/Units::C << std::endl;
  std::cout << "inclination = " << incl << " degrees" << std::endl;
  
  std::cout << std::endl;
 
  // Units conversions.
  incl *= Units::PI / 180.0;
  theta *= Units::PI / 180.0;
  mu = cos(theta);
  mass = Units::cgs_to_nounits( mass*Units::MSUN, Units::MASS );
  req = Units::cgs_to_nounits( req*1.0e5, Units::LENGTH );
  omega = Units::cgs_to_nounits( 2.0*Units::PI*omega, Units::INVTIME );
  distance = Units::cgs_to_nounits( distance*100, Units::LENGTH );

  double baromega; // dimensionless omega
  baromega = omega * sqrt( pow(req,3)/mass);
  std::cout << "baromega = " << baromega <<std::endl;

  north = incl;
  south = Units::PI - incl;

  std::cout << "mass/radius = " << mass/req << std::endl;

  // curve is a structure holding parameters describing the star and emission properties
  curve.para.mass = mass;
  curve.para.omega = omega;
  //  curve.para.theta = theta;
  //curve.para.incl = incl;
  curve.para.aniso = aniso;
  curve.para.Gamma = Gamma;
  curve.para.bbrat = bbrat;
  curve.para.ts = ts;

 
  // Set up model describing the oblate shape of the star
  OblModelBase* model;
  if (modelchoice == 1 ) {
    // Default model 
    model = new PolyOblModelNHQS( rspot, req,
				  PolyOblModelBase::zetaparam(mass,req),
				  PolyOblModelBase::epsparam(omega, mass, req)
				  );
  }
  else if(modelchoice == 2 ) {
    // Alternative model for quark stars (not very different)
    model = new PolyOblModelCFLQS( rspot, req,
				   PolyOblModelBase::zetaparam(mass,req),
				   PolyOblModelBase::epsparam(omega, mass, req)
				   );
  }
  else if(modelchoice == 3 ) {
    // Use a spherical star
    model = new SphericalOblModel( rspot );
  }
  else {
    std::cerr << "Invalid modelchoice parameter. Exiting." << std::endl;
    return -1;
  }

  // defltoa is a structure that "points" to routines in the file "OblDeflectionTOA.cpp"
  // used to compute deflection angles and times of arrivals 
  OblDeflectionTOA* defltoa = new OblDeflectionTOA( model, mass );

 

  curve.junk.shift_t = shift_t;
  curve.junk.infile_is_set = infile_is_set;
  curve.numbins = numbins;

  //  std::cout << "i" << "\t" << "a1" << "\t" << "phase" << "\t" << "Chi-Square" << std::endl;

  double dOmega_mu(0.0), dFlux_mu(0.0);
  double Omega_s(0.0), Flux(0.0);
  double Omega_sQ(0.0), Flux_Q(0.0);
  double dBoloFlux_mu(0.0), BoloFlux(0.0);
  double dArea_mu(0.0), Area(0.0);

  double red; // Gravitational redshift red=1+z= (1-2M/R)^{-1/2}
  double quad; // Quadrupole Moment
  double riso; // Isotropic radius
  double framedragging; 
  double inertia;
  double speed;

 

  if (quadmodel == 1){
    // quad = - 0.1 * pow(baromega,2) * pow( req/mass, 2);
    quad = -0.2;
  }
  else 
    quad = 0.0;
  std::cout << "mass/req = " << mass/req 
    << "Quadrupole = " << quad  << std::endl;

  double a_kerr;

  if (fdmodel==1)
    a_kerr=0.24;
  else
    a_kerr=0.0;


  if (fdmodel==1)
    inertia = sqrt(mass/req) * (1.14 - 2.53*mass/req + 5.6*pow(mass/req,2));
  else
    inertia = 0;

  double dmu,dphi, P2;
  dmu = 1.0/(nmu*1.0);
  dphi = 2 * Units::PI / (numbins*1.0);
 

  // Loop through the star's latitudes
   for (int j(0); j<nmu; j++){

     mu = j*dmu;
     theta = acos(mu);
     curve.para.theta = theta;
     //     P2 = 0.5 * ( 3.0*mu*mu - 1.0);

     // Values we need in some of the formulas.

     if (modelchoice==1){
       rspot = calcrspot( omega, mass, theta, req);
       cosg = cosgamma(mu,req,rspot,mass,omega);
     }
     if (modelchoice==3){
       rspot = req;
       cosg = 1.0;
     }

     //     curve.para.dS = pow(rspot/distance,2)/cosg * dmu * dphi;  

     curve.para.dS = dmu * dphi/cosg;     
     curve.para.cosgamma = cosg;
     curve.para.radius = rspot;

     riso = rspot * pow( 1 + 0.5*mass/rspot, -2);
     //     red = 1.0/sqrt(1 - 2*mass/rspot) * exp( quad * P2 * pow( mass/riso,3)) ;

     double x;
     x = rspot/mass;

     double F1;
     F1 = -5.0/8.0*(x-1.0)/(x*(x-2)) * (2.0+6.0*x - 3.0*x*x)
       - 15.0/16.0 * x*(x-2) * log(x/(x-2.0));

     double Sigma;
     Sigma = pow(x,2) + pow(a_kerr*mu,2);

     double Nsquare;
     Nsquare = (1.0 - 2.0*x/Sigma) * ( 1.0 + 2.0*quad*F1*P2);

     red = 1.0/sqrt(Nsquare);

     std::cout
       << "j = " << j
       << " mu = " << mu 
       << " dS = " << curve.para.dS
       <<  std::endl;

     framedragging = 2.0*inertia * pow(req/mass,2) * pow(mass/riso,3) * (1.0 - 3.0*mass/riso);
     speed = omega*rspot*sin(theta)* red * (1.0 - framedragging);

     //     std::cout << "mu = " << mu << " P2(mu) = " << P2 << " red^{-4} = " << pow(red*sqrt(1 - 2/x),-4) << std::endl;
     //std::cout << "mass/radius = " << mass/req << " Quad = " << quad << std::endl;
     //std::cout << "inertia = " << inertia << " fd = " << framedragging << " speed = " << speed << std::endl;


    
       
     curve.defl.b_max =  defltoa->bmax_outgoing(rspot);
     curve.defl.psi_max = defltoa->psi_max_outgoing(curve.defl.b_max,rspot,&curve.problem);

       // Compute b vs psi lookup table good for the specified M/R and mu
     b_mid = curve.defl.b_max*0.9;
     curve.defl.b_psi[0] = 0.0;
     curve.defl.psi_b[0] = 0.0;
     for (unsigned int i(1); i < NN+1; i++){ // compute table of b vs psi points 
       curve.defl.b_psi[i] = b_mid * i/(NN*1.0);
       curve.defl.psi_b[i] = defltoa->psi_outgoing(curve.defl.b_psi[i],rspot,curve.defl.b_max, curve.defl.psi_max, &curve.problem);
     }
     // For arcane reasons, the table is not evenly spaced.
     for (unsigned int i(NN+1); i < 3*NN; i++){ // compute table of b vs psi points 
       curve.defl.b_psi[i] = b_mid + (curve.defl.b_max-b_mid)/2.0 * (i-NN)/(NN*1.0);
       curve.defl.psi_b[i] = defltoa->psi_outgoing(curve.defl.b_psi[i],rspot,curve.defl.b_max, curve.defl.psi_max, &curve.problem);
     }
     curve.defl.b_psi[3*NN] = curve.defl.b_max;
     curve.defl.psi_b[3*NN] = curve.defl.psi_max;
     // Finished computing lookup table
     
     //std::cout << "Finished computing lookup table" << std::endl;


     // Loop through the hemispheres
     for (int k(0);k<=1; k++){

       if (k==1) curve.para.incl = south;
       else curve.para.incl = north;

       curve =  ComputeAngles( &curve, defltoa);

       dOmega_mu = 0.0;
       dBoloFlux_mu = 0.0;
       dFlux_mu = 0.0;
       dArea_mu = 0.0;

       double eta;
       double redshift;

       redshift = 1.0/sqrt(1 - 2/x);

       for (int i(0);i<numbins;i++){

	 dArea_mu += curve.para.dS * pow( Units::nounits_to_cgs( rspot*1.0e-5, Units::LENGTH ),2) ;

	 if (curve.dOmega_s[i] != 0.0){

	   eta = sqrt( 1 - speed*speed)/(1.0 - speed * curve.cosxi[i]);
	   
	   dOmega_mu += curve.dOmega_s[i];

	   dBoloFlux_mu += curve.dOmega_s[i]*pow(eta,4)/Units::PI;

	   dFlux_mu += curve.dOmega_s[i] * pow(eta,3) * BlackBody(Temperature,E0*redshift/curve.eta[i]);
	   
	 }
       }

       //std::cout << "dA = " << dArea_mu << std::endl;
       //std::cout << "dFlux_mu = " << dFlux_mu << std::endl;

       
       Area += dArea_mu  ;

       Omega_s += dOmega_mu * pow(rspot/distance,2);
       Flux += dFlux_mu * pow( sqrt(1 - 2/x), 3) * (1.0 / ( E0 * Units::H_PLANCK )) * pow(rspot/distance,2);
       BoloFlux += dBoloFlux_mu * pow( red,-4) * pow(rspot/distance,2);

    

     }
   }

   std::cout << "Area = " << Area << std::endl;

  std::cout << "pi R^2/D^2 (1+z)^2 = " <<   Units::PI * pow(req,2)/(1.0-2.0*mass/req) << std::endl;

  std::cout << "R^2/D^2 (1+z)^{-2} = " << pow(req,2)*(1.0-2.0*mass/req) << std::endl;

 // Print out information about the model
  std::cout << std::endl << "Pulse: Stellar Model Parameters" << std::endl << std::endl;
  std::cout << "Mass = " << Units::nounits_to_cgs( mass/Units::MSUN, Units::MASS ) << " Msun" << std::endl;
  std::cout << "Req = " << Units::nounits_to_cgs( req*1.0e-5, Units::LENGTH ) << " km" <<  std::endl;
  std::cout << "2GM/Rc^2 = " << 2.0*Units::G*mass*Units::MSUN/(rspot*1e+5*Units::C*Units::C)  <<  std::endl;
  std::cout << "spin = " << omega_cgs << " Hz" << std::endl;
  std::cout << "v/c = " << rspot * 1e6 * omega_cgs * 2.0*Units::PI/Units::C << std::endl;
  std::cout << "inclination = " << incl << " degrees" << std::endl;
  
  std::cout << std::endl;

  out.open(out_file, std::ios::app);

  //  out = fopen(out_file, "a");
	
 // Print out information about the model
  out << "#M = " << Units::nounits_to_cgs( mass/Units::MSUN, Units::MASS ) << " Msun" 
      << " #Req = " << Units::nounits_to_cgs( req*1.0e-5, Units::LENGTH ) << " km" << std::endl;
    
	  
  out << "#spin      "
      << "Flux(1keV) "  
      << "BolFlux/Fs "
      << "v_{eq}/c   "
      << "Area       "
      << "A/4piR^2   "
      << "SolidAng   "
      << "SolidAng/Ss" << std::endl;

  out  << omega_cgs
       << "          " <<  Flux 
       << "          " <<  BoloFlux/(pow(req/distance,2)*(1.0-2.0*mass/req)) 
       << "          " <<     req * 1e6 * omega_cgs * 2.0*Units::PI/Units::C * red
       << "          " << Area
       << "          " << Area/(4.0 * Units::PI * pow( Units::nounits_to_cgs( rspot*1.0e-5, Units::LENGTH ),2))
       << "          " << Omega_s
       << "          " << Omega_s/(Units::PI * pow(req/distance,2)/(1.0-2.0*mass/req))
      << std::endl;



  delete defltoa;
  delete model;
  return 0;

 } catch(std::exception& e) {
  std::cerr << "Top-level exception caught in application oblflux: " << std::endl
	    << e.what() << std::endl;
  return -1;
 }
