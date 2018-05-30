// SphericalOblModel.h
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

#ifndef SPHERICALOBLMODEL_H
#define SPHERICALOBLMODEL_H

#include "OblModelBase.h"
#include <exception>

class SphericalOblModel : public OblModelBase {
 private:
  const double Req;
 public:
  SphericalOblModel( const double& Req_val );
  double R_at_costheta( const double& costheta ) const throw(std::exception);
  double Dtheta_R( const double& costheta ) const throw(std::exception);
  double f(const double& costheta) const throw(std::exception);
  double cos_gamma(const double& costheta) const throw(std::exception);
};

#endif // SPHERICALOBLMODEL_H
