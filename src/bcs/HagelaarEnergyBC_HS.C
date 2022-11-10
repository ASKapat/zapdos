//* This file is part of Zapdos, an open-source
//* application for the simulation of plasmas
//* https://github.com/shannon-lab/zapdos
//*
//* Zapdos is powered by the MOOSE Framework
//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HagelaarEnergyBC_HS.h"

registerADMooseObject("ZapdosApp", HagelaarEnergyBC_HS);

InputParameters
HagelaarEnergyBC_HS::validParams()
{
  InputParameters params = ADIntegratedBC::validParams();
  params.addRequiredParam<Real>("r", "The reflection coefficient");
  params.addRequiredCoupledVar("potential", "The electric potential");
  params.addRequiredCoupledVar("HS_density", "The heavy species density.");
  params.addRequiredParam<Real>("position_units", "Units of position.");
  return params;
}

HagelaarEnergyBC_HS::HagelaarEnergyBC_HS(const InputParameters & parameters)
  : ADIntegratedBC(parameters),
    _r_units(1. / getParam<Real>("position_units")),
    _r(getParam<Real>("r")),

    // Coupled Variables
    _grad_potential(adCoupledGradient("potential")),
    _n_HS(adCoupledValue("HS_density")),
    
    _muip(getADMaterialProperty<Real>("mu" + (*getVar("HS_density", 0)).name())),

    _massip(getMaterialProperty<Real>("mass" + (*getVar("HS_density", 0)).name())),
    _sgnip (getMaterialProperty<Real>("sgn" + (*getVar("HS_density", 0)).name())),

    _e(getMaterialProperty<Real>("e")),

    _muE_HS(getADMaterialProperty<Real>("mu" + _var.name()))

{
  _a = 0.5;
  _v_thermal = 0.0;
}

ADReal
HagelaarEnergyBC_HS::computeQpResidual()
{
  if (_normals[_qp] * _sgnip[_qp]* -_grad_potential[_qp] > 0.0)
  {
    _a = 1.0;
  }
  else
  {
    _a = 0.0;
  }
  _v_thermal =
      std::sqrt(8 * _e[_qp] * 2.0 / 3 * std::exp(_u[_qp] - _n_HS[_qp]) / (M_PI * _massip[_qp]));

  return _test[_i][_qp] * _r_units * (1. - _r) / (1. + _r) *
         (_sgnip[_qp]*(2. * _a - 1.) * _muE_HS[_qp] * -_grad_potential[_qp] * _r_units * _normals[_qp] +
          5. / 6. * _v_thermal) *
         std::exp(_u[_qp]);
}
