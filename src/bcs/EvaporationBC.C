#include "EvaporationBC.h"

registerADMooseObject("ZapdosApp", EvaporationBC);

InputParameters
EvaporationBC::validParams()
{
InputParameters params = ADIntegratedBC::validParams();
  params.addRequiredCoupledVar("Wall_temperature","Temperature of wall");
  params.addRequiredParam<Real>("position_units", "Units of position.");
  params.addRequiredParam<Real>("Vapor_pressure_const", "Constant term in polynomial vapor pressure equation");
  params.addRequiredParam<Real>("Vapor_pressure_inverse_coeff", "Inverse T coefficient in polynomial vapor pressure equation");
  params.addClassDescription("Evaporative flux based on Langmuir evaporation law,"
                             "Vapor pressure equations from Alcock et al.,"
                             "Canadian Metallurgical Quarterly 23 (1984)");


  return params;
}

EvaporationBC::EvaporationBC(const InputParameters & parameters)
	: ADIntegratedBC(parameters),
	  _T(adCoupledValue("Wall_temperature")),
      _r_units(1. / getParam<Real>("position_units")),
      _a(getParam<Real>("Vapor_pressure_const")),
      _b(getParam<Real>("Vapor_pressure_inverse_coeff")),
      _m(getMaterialProperty<Real>("mass" + _var.name())),
      _kb(getMaterialProperty<Real>("k_boltz"))
{
}
ADReal
EvaporationBC::computeQpResidual()
{
_evap_flux = std::pow(10, _a-_b/_T[_qp])/(std::sqrt(2*M_PI*_m[_qp]*_kb[_qp]*_T[_qp]));


  return -_test[_i][_qp] * _r_units  * _evap_flux;
}