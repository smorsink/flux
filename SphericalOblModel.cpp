// SphericalOblModel.cpp
//
// For Spherical Stars
// (C) Coire Cadeau, 2007

// Source (C) Coire Cadeau 2007, all rights reserved.
//
// Permission is granted for private use only, and not
// distribution, either verbatim or of derivative works,
// in whole or in part.
//
// The code is not thoroughly tested or guaranteed for
// any particular use.

#include "SphericalOblModel.h"

SphericalOblModel::SphericalOblModel( const double& Req_val ) : Req(Req_val) { }

double SphericalOblModel::R_at_costheta( const double& costheta ) const throw(std::exception) { 
  return Req; 
}

double SphericalOblModel::Dtheta_R( const double& costheta ) const throw(std::exception) { 
  return double(0.0);
}

double SphericalOblModel::f(const double& costheta) const throw(std::exception) { 
  return double(0.0); 
}

double SphericalOblModel::cos_gamma(const double& costheta) const throw(std::exception) { 
  return double(1.0); 
}
