//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HagelaarIonDiffusionBC_Ti_Te.h"

registerADMooseObject("ZapdosApp", HagelaarIonDiffusionBC_Ti_Te);

InputParameters
HagelaarIonDiffusionBC_Ti_Te::validParams()
{
  InputParameters params = ADIntegratedBC::validParams();
  params.addRequiredParam<Real>("r", "The reflection coefficient");
  params.addRequiredParam<Real>("position_units", "Units of position.");
  params.addRequiredCoupledVar("mean_en", "The mean energy.");
  params.addParam<Real>(
      "user_velocity", -1., "Optional parameter if user wants to specify the thermal velocity.");
  params.addClassDescription("Kinetic electron boundary condition"
                             "(Based on DOI:https://doi.org/10.1103/PhysRevE.62.1452)");
  return params;
}

HagelaarIonDiffusionBC_Ti_Te::HagelaarIonDiffusionBC_Ti_Te(const InputParameters & parameters)
  : ADIntegratedBC(parameters),
    _r_units(1. / getParam<Real>("position_units")),
    _r(getParam<Real>("r")),
    
    // Coupled Variables
    _mean_en(adCoupledValue("mean_en")),

    _kb(getMaterialProperty<Real>("k_boltz")),
    _T(getADMaterialProperty<Real>("T" + _var.name())),
    _mass(getMaterialProperty<Real>("mass" + _var.name())),
    _e(getMaterialProperty<Real>("e")),
    _user_velocity(getParam<Real>("user_velocity"))
{
  _v_thermal = 0.0;
}

ADReal
HagelaarIonDiffusionBC_Ti_Te::computeQpResidual()
{
  if (_user_velocity > 0.)
    _v_thermal = _user_velocity;
  else
  	_v_thermal =
      std::sqrt(8 * _e[_qp] * 2.0 / 3 * std::exp(_mean_en[_qp] - _u[_qp]) / (M_PI * _mass[_qp]));
    //*_v_thermal = std::sqrt(8 * _kb[_qp] * _T[_qp] / (M_PI * _mass[_qp]));

  return _test[_i][_qp] * _r_units * (1. - _r) / (1. + _r) * 0.5 * _v_thermal * std::exp(_u[_qp]);
}
