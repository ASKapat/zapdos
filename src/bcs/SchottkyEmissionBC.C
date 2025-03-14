//* This file is part of Zapdos, an open-source
//* application for the simulation of plasmas
//* https://github.com/shannon-lab/zapdos
//*
//* Zapdos is powered by the MOOSE Framework
//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SchottkyEmissionBC.h"

registerMooseObject("ZapdosApp", SchottkyEmissionBC);

InputParameters
SchottkyEmissionBC::validParams()
{
  InputParameters params = ADIntegratedBC::validParams();
  params.addRequiredParam<Real>("r", "The reflection coefficient");
  params.addRequiredCoupledVar("potential", "The electric potential");
  params.addRequiredCoupledVar("mean_en", "The mean energy.");
  params.addRequiredCoupledVar("ip", "The ion density.");
  params.addRequiredParam<Real>("position_units", "Units of position.");
  params.addRequiredParam<std::string>("potential_units", "The potential units.");
  params.addParam<Real>("tau", 1e-9, "The time constant for ramping the boundary condition.");
  params.addParam<bool>("relax", false, "Use relaxation for emission.");
  return params;
}

SchottkyEmissionBC::SchottkyEmissionBC(const InputParameters & parameters)
  : ADIntegratedBC(parameters),

    _r_units(1. / getParam<Real>("position_units")),
    _r(getParam<Real>("r")),

    // Coupled Variables
    _grad_potential(adCoupledGradient("potential")),
    _mean_en(adCoupledValue("mean_en")),
    _ip_var(*getVar("ip", 0)),
    _ip(adCoupledValue("ip")),
    _grad_ip(adCoupledGradient("ip")),

    _massem(getMaterialProperty<Real>("massem")),
    _e(getMaterialProperty<Real>("e")),
    _sgnip(getMaterialProperty<Real>("sgn" + _ip_var.name())),
    _muip(getADMaterialProperty<Real>("mu" + _ip_var.name())),
    _Dip(getADMaterialProperty<Real>("diff" + _ip_var.name())),
    _se_coeff(getMaterialProperty<Real>("se_coeff")),
    _work_function(getMaterialProperty<Real>("work_function")),
    _field_enhancement(getMaterialProperty<Real>("field_enhancement")),
    _Richardson_coefficient(getMaterialProperty<Real>("Richardson_coefficient")),
    _cathode_temperature(getMaterialProperty<Real>("cathode_temperature")),
    _a(0.5),
    _v_thermal(0),
    _ion_flux(0, 0, 0),
    _tau(getParam<Real>("tau")),
    _relax(getParam<bool>("relax")),
    _potential_units(getParam<std::string>("potential_units"))

{
  if (_potential_units.compare("V") == 0)
  {
    _voltage_scaling = 1.;
  }
  else if (_potential_units.compare("kV") == 0)
  {
    _voltage_scaling = 1000;
  }

  _dPhi_over_F =
      3.7946865E-5 * std::sqrt(_voltage_scaling); // eV*sqrt(m/kV) (if _voltage_scaling = 1000)
}

ADReal
SchottkyEmissionBC::computeQpResidual()
{
  ADReal dPhi;
  Real kB;
  ADReal jRD;
  ADReal jSE;
  ADReal F;
  Real _relaxation_Expr;

  _v_thermal =
      std::sqrt(8 * _e[_qp] * 2.0 / 3 * std::exp(_mean_en[_qp] - _u[_qp]) / (M_PI * _massem[_qp]));

  if (_normals[_qp] * -1.0 * -_grad_potential[_qp] > 0.0)
  {
    _a = 1.0;
    return 0;
  }
  else
  {
    _a = 0.0;

    _ion_flux = _sgnip[_qp] * _muip[_qp] * -_grad_potential[_qp] * _r_units * std::exp(_ip[_qp]) -
                _Dip[_qp] * std::exp(_ip[_qp]) * _grad_ip[_qp] * _r_units;

    // Schottky emission
    // je = AR * T^2 * exp(-(wf - dPhi) / (kB * T))
    // dPhi = _dPhi_over_F * sqrt(F) // eV

    F = -(1 - _a) * _field_enhancement[_qp] * _normals[_qp] * _grad_potential[_qp] * _r_units;

    kB = 8.617385E-5; // eV/K
    dPhi = _dPhi_over_F * std::sqrt(F);

    jRD = _Richardson_coefficient[_qp] * std::pow(_cathode_temperature[_qp], 2) *
          std::exp(-(_work_function[_qp] - dPhi) / (kB * _cathode_temperature[_qp]));
    jSE = _e[_qp] * 6.02E23 * _se_coeff[_qp] * _ion_flux * _normals[_qp];

    if (_relax)
    {
      _relaxation_Expr = std::tanh(_t / _tau);
    }
    else
    {
      _relaxation_Expr = 1.0;
    }

    return _relaxation_Expr * _test[_i][_qp] * _r_units * 2. / (1. + _r) * (1 - _a) * (-jRD - jSE) /
           (_e[_qp] * 6.02E23);
  }
}
