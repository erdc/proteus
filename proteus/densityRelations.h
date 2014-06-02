#ifndef DENSITYRELATIONS_H
#define DENSITYRELATIONS_H

/** \file densityRelations.h
    \defgroup densityrelations densityrelations
    \brief A library of density-pressure relations
    @{
*/

/** The density of  an incompressible fluid.

    @param[in] psi The pressure head (not used)
    @param[in] rwork[0] The density constant
    @param[out] rho The density
    @param[out] drho The derivative of density with respect to pressure head

    This is a trivial function but it demonstrates the density relation interface.
    \f[ \rho(\psi) = \rho_0 \f]
*/
class DensityRelation
{
public:
  double rho,drho;
};

class ConstantDensity : public DensityRelation
{
public:
  double rho_0;
  inline ConstantDensity(const double* rwork):
    rho_0(rwork[0])
  {
    rho=rho_0;
    drho=0.0;
  }
  inline void calc(const double& psi)
  {
    rho=rho_0;
    drho=0.0;
  }
};

/** The density of a linearly compressible fluid.

    @param[in] psi The pressure head (not used)
    @param[in] rwork[0] The reference density, \f$ \rho_0 \f$
    @param[in] rwork[1] The reference pressure head, \f$ \psi_0 \f$
    @param[in] rwork[2] The compressibility \f$ \beta \f$
    @param[out] rho The density
    @param[out] drho The derivative of density with respect to pressure head

    This is a trivial function but it demonstrates the density relation interface.
    \f[ \rho &= \rho_0 [1+\beta(\psi - \psi_0) ] \f]
*/
class LinearDensity : public DensityRelation
{
public:
  double rho_0,psi_0,beta;
  inline LinearDensity(const double* rwork):
    rho_0(rwork[0]),
    psi_0(rwork[1]),
    beta(rwork[2])
  {}
  
  inline void calc(const double& psi)
  {
    rho = rho_0*(1.0+beta*(psi - psi_0));
    drho = beta*rho_0;
  }
};

/** The density of an exponentially compressible fluid.

    @param[in] psi The pressure head (not used)
    @param[in] rwork[0] The reference density, \f$ \rho_0 \f$
    @param[in] rwork[1] The reference pressure head, \f$ \psi_0 \f$
    @param[in] rwork[2] The compressibility \f$ \beta \f$
    @param[out] rho The density
    @param[out] drho The derivative of density with respect to pressure head

    This is a trivial function but it demonstrates the density relation interface.
    \f[ \rho &= \rho_0 e^{\beta_f \psi} \f]
*/
class ExponentialDensity : public LinearDensity
{
public:
  inline ExponentialDensity(const double* rwork):
    LinearDensity(rwork)
  {}
  inline void calc(const double& psi)
  {
    rho = rho_0*exp(beta*(psi-psi_0));
    drho = beta*rho;
  }
};

/** The density of a real gas.

    @param[in] psi The pressure head (not used)
    @param[in] rwork[0] The reference pressure head, \f$ \psi_0 \f$
    @param[in] rwork[1] The constant, \f$ 1/ZRT \f$, where \f$ Z \f$ is the real gas constant, \f$ R \f$ is the ideal gas constant, and \f$ T \f$ is the temperature.
    @param[out] rho The density
    @param[out] drho The derivative of density with respect to pressure head

    \f[ \rho &= (\psi - \psi_0)/ZRT \f]
*/
class RealGasDensity : public DensityRelation
{
public:
  double psi_0,oneOverZRT;
  inline RealGasDensity(const double *rwork):
    psi_0(rwork[0]),
    oneOverZRT(rwork[1])
  {}
  inline void calc(const double& psi)
  {
    rho = psi*oneOverZRT;
    drho = oneOverZRT;
  }
};

/** The density of an ideal gas.

    @param[in] psi The pressure head
    @param[in] rwork[0] The temperature in say [K], \f$ T \f$
    @param[in] rwork[1] The molar weight of the gas in say [kg]/[mol], \f$ W \f$,
    @param[in] rwork[2] The ideal gas weight constant in say [J]/[mol]/[K], \f$ R \f$
    @param[in] rwork[3] Conversion factor for input pressure head to pressure units consistent with W, R, and T, ie rho_0*|g|
    @param[in] rwork[4] reference density, say in [kg]/[m^3]
    @param[in] rwork[5] reference pressure head, say in [m]
    @param[out] rho The density
    @param[out] drho The derivative of density with respect to pressure head

    \f[ \rho &= rho_0 + (\psi - psi0) d W/RT \f]
*/
class IdealGasDensity : public DensityRelation
{
public:
  double T,W,R,convFactor,WoverRT,rho0,psi0;
  inline IdealGasDensity(const double *rwork):
    T(rwork[0]),
    W(rwork[1]),
    R(rwork[2]),
    convFactor(rwork[3]),
    rho0(rwork[4]),
    psi0(rwork[5])
    {
      WoverRT = convFactor*W/(R*T);
    }
  inline void calc(const double& psi)
  {
    rho = rho0 + (psi-psi0)*WoverRT;
    drho = WoverRT;
  }
};


/** @} */
#endif
