//* This file is part of Zapdos, an open-source
//* application for the simulation of plasmas
//* https://github.com/shannon-lab/zapdos
//*
//* Zapdos is powered by the MOOSE Framework
//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FloatingPotentialNeumann.h"

// MOOSE includes
#include "Function.h"

registerADMooseObject("ZapdosApp", FloatingPotentialNeumann);

InputParameters
FloatingPotentialNeumann::validParams()
{
  InputParameters params = ADIntegratedBC::validParams();
  params.addRequiredCoupledVar("em", "Electron density");
  params.addRequiredCoupledVar("mean_en", "Electron energy");
  params.addRequiredCoupledVar("ip", "Ion species(es) incident on surface");
  params.addRequiredCoupledVar("ion_energies", "Energy of ion(s) upstream");
  params.addRequiredParam<Real>("position_units", "Units of position");  
  params.addRequiredParam<std::string>("potential_units", "The potential units.");
  return params;
}

FloatingPotentialNeumann::FloatingPotentialNeumann(const InputParameters & parameters)
  : ADIntegratedBC(parameters),

    _em(adCoupledValue("em")),
    _mean_en(adCoupledValue("mean_en")),
    _potential_units(getParam<std::string>("potential_units")),
    _r_units(1. / getParam<Real>("position_units")),
    _eps(getMaterialProperty<Real>("eps")),
    _N_A(getMaterialProperty<Real>("N_A")),
    _muem(getADMaterialProperty<Real>("muem")),
    _e(getMaterialProperty<Real>("e")),
    _massem(getMaterialProperty<Real>("massem")),
    _kb(getMaterialProperty<Real>("k_boltz"))
{

  _a_e = 0.0;
  
  if (_potential_units.compare("V") == 0)
    _voltage_scaling = 1.;
  else if (_potential_units.compare("kV") == 0)
    _voltage_scaling = 1000;

  _num_ions = coupledComponents("ip");

  // Resize the vectors to store _num_ions values:
  _ip.resize(_num_ions);
  _Eip.resize(_num_ions);
  _muip.resize(_num_ions);
  _sgnip.resize(_num_ions);
  _massip.resize(_num_ions);
  _a.resize(_num_ions);

  // Retrieve the values for each ion and store in the relevant vectors.
  // Note that these need to be dereferenced to get the values inside the
  // main body of the code.
  // e.g. instead of "_ip[_qp]" it would be "(*_ip[i])[_qp]", where "i"
  // refers to a single ion species.
  for (unsigned int j = 0; j < _ip.size(); ++j)
  {
    _ip[j] = &adCoupledValue("ip", j);
    _Eip[j] = &adCoupledValue("ion_energies", j);
    _muip[j] = &getADMaterialProperty<Real>("mu" + (*getVar("ip", j)).name());
    _sgnip[j] = &getMaterialProperty<Real>("sgn" + (*getVar("ip", j)).name());
    _massip[j] = &getMaterialProperty<Real>("mass" + (*getVar("ip", j)).name());
    _a[j] = 0.0;
    
  }
}

ADReal
FloatingPotentialNeumann::computeQpResidual()
{

  if (_normals[_qp] * -1.0 * -_grad_u[_qp] > 0.0)
  {
    _a_e = 0.0;
  }
  else
  {
    _a_e = 1.0;
  }

  _ion_therm_drift = 0.0;
  _advect = 0.0;

  for (unsigned int k = 0; k < _num_ions; ++k)
  {
    if (_normals[_qp] * (*_sgnip[k])[_qp] * -_grad_u[_qp] > 0.0)
    {
    _a[k] = 1.0;
  }
    else
    {
    _a[k] = 0.0;
  }
    _v_th_i = std::sqrt(8 * _e[_qp] * 2.0 / 3 * std::exp((*_Eip[k])[_qp] - (*_ip[k])[_qp])
        / (M_PI * (*_massip[k])[_qp]));
        
    _ion_therm_drift += std::exp((*_ip[k])[_qp])*_v_th_i;
    
    _advect += (2*_a[k]-1)*std::abs((*_muip[k])[_qp])*std::exp((*_ip[k])[_qp]);
  }
  _v_e_th = std::sqrt(8 * _e[_qp] * 2.0 / 3 * std::exp(_mean_en[_qp] - _em[_qp]) /
                      (M_PI * _massem[_qp]));
                      
//   std::cout <<(_ion_therm_drift - _v_e_th*std::exp(_em[_qp]))/
//                       (_advect+(2*_a_e-1)*_muem[_qp]*std::exp(_em[_qp]))<< std::endl;
                      
  return -_test[_i][_qp] * _eps[_qp] * (_ion_therm_drift - _v_e_th*std::exp(_em[_qp]))/
                      (_advect+(2*_a_e-1)*_muem[_qp]*std::exp(_em[_qp]));
}