//* This file is part of Zapdos, an open-source
//* application for the simulation of plasmas
//* https://github.com/shannon-lab/zapdos
//*
//* Zapdos is powered by the MOOSE Framework
//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FloatingPotential.h"
#include "Function.h"

registerMooseObject("ZapdosApp", FloatingPotential);

InputParameters
FloatingPotential::validParams()
{
  InputParameters params = ADIntegratedBC::validParams();
  
  params.addRequiredParam<Real>("V_plasma", "Plasma potential");
  
  params.addRequiredCoupledVar("em", "Electron (density)");
  params.addRequiredCoupledVar("mean_en", "Energy of electrons");
  params.addRequiredParam<std::vector<Real>>("mass_ion", "Ion masses");
  params.addRequiredCoupledVar("ion", "Ion species(es) incident on surface");
  params.addRequiredCoupledVar("ion_energies", "Energy of ion(s) upstream");
  params.addRequiredParam<Real>("position_units", "Units of position");  
  params.addRequiredParam<std::string>("potential_units", "The potential units.");
  params.addRequiredParam<Real>("penalty", "Penalty value");

  params.addClassDescription("Dirichlet penalty boundary condition for floating potential");
  return params;
}

FloatingPotential::FloatingPotential(const InputParameters & parameters)
  : ADIntegratedBC(parameters),

    _Vp(getParam<Real>("V_plasma")),
    _em(coupledValue("em")),
    _Te(coupledValue("mean_en")),
    _massip(getParam<std::vector<Real>>("mass_ion")),
    _potential_units(getParam<std::string>("potential_units")),
    _r_units(1. / getParam<Real>("position_units"))

{
    if (_potential_units.compare("V") == 0)
    _voltage_scaling = 1.;
    else if (_potential_units.compare("kV") == 0)
     _voltage_scaling = 1000;

  _num_ions = coupledComponents("ion");
  _n_ip.resize(_num_ions);
  _Eip.resize(_num_ions);

  for (unsigned int i = 0; i < _num_ions; ++i)
  {
    _n_ip[i] = &coupledValue("ion", i);
    _Eip[i] = &coupledValue("ion_energies", i);
  }
  
  
}

ADReal
FloatingPotential::computeQpResidual()
{

// _ion_flux = 0
// _ion_flux_tot = 0
// for (unsigned int i = 0; i < _num_ions; ++i)
//   { 
//   _ion_flux =   std::exp(_n_ip[i])*std::sqrt(std::exp((*_Eip[i])[_qp] - (*_n_ip[i])[_qp])
//      / (*_massip[i])[_qp])
//   _ion_flux_tot += _ion_flux
//   }
// return std::exp((_mean_en[_qp] - _em[_qp])*(std::log(std::sqrt(2*M_PI*_massem[_qp]/(std::exp(_Te[_qp] - _em[_qp])))/_em[_qp])/_voltage_scaling
//   +std::log(_ion_flux_tot))+_Vp;
  
_ion_flux = 0;
_ion_flux_tot = 0;
for (unsigned int i = 0; i < _num_ions; ++i)
//   {
//   _ion_flux =   std::exp((*_n_ip[i])[4068])*std::sqrt(std::exp((*_Eip[i])[4068] - (*_n_ip[i])[4068])/_massip[i]);
//   _ion_flux_tot += _ion_flux;
//   }
// _V_f =  std::exp(_Te[4068] - _em[4068])*(std::log(std::sqrt(2*M_PI*(9.11E-31)/(std::exp(_Te[4068] - _em[4068])))/_em[4068])
//   +std::log(_ion_flux_tot))/_voltage_scaling + _Vp;
//   
//   return _p * _test[_i][_qp] * (-_V_f + _u[_qp]);

  {
//   _ion_flux = std::exp(_n_ip[i])*std::sqrt(std::exp(_Eip[i] - _n_ip[i])/_massip[i]);
//   _ion_flux_tot += _ion_flux;
  _ion_flux =   std::exp((*_n_ip[i])[_qp])*std::sqrt(std::exp((*_Eip[i])[_qp] - (*_n_ip[i])[_qp])/_massip[i]);
  _ion_flux_tot += _ion_flux;
  }

_V_f =  std::exp(_Te[_qp] - _em[_qp])*(std::log(std::sqrt(2*M_PI*(9.11E-31)/(std::exp(_Te[_qp] - _em[_qp])))/_em[_qp])
  +std::log(_ion_flux_tot))/_voltage_scaling + _Vp;
  
  return _p * _test[_i][_qp] * (-_V_f + _u[_qp]);
}