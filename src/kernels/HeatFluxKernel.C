//* This file is part of Zapdos, an open-source
//* application for the simulation of plasmas
//* https://github.com/shannon-lab/zapdos
//*
//* Zapdos is powered by the MOOSE Framework
//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeatFluxKernel.h"

registerADMooseObject("ZapdosApp", HeatFluxKernel);

InputParameters
HeatFluxKernel::validParams()
{
  InputParameters params = ADKernel::validParams();

  params.addRequiredCoupledVar("potential", "The potential that drives the advective flux.");
  params.addRequiredParam<Real>("position_units", "Units of position.");
  params.addRequiredCoupledVar("incident", "All of the  species incident with this boundary.");
  params.addRequiredCoupledVar("incident_energies",
                               "Corresponding energies of species incident with this boundary.");

  params.addClassDescription("Calculates heat flux for a Variable");
  return params;
}

HeatFluxKernel::HeatFluxKernel(const InputParameters & parameters)
  : ADKernel(parameters),
  
    _r_units(1. / getParam<Real>("position_units")),
    _grad_potential(adCoupledGradient("potential"))
    
    
{

  // Resize the vectors to store _num_species pointers
  _num_inc = coupledComponents("incident");
  _n.resize(_num_inc);
  _mu.resize(_num_inc);
  _muE.resize(_num_inc);
  _E.resize(_num_inc);
  _sgn.resize(_num_inc);
  _mass.resize(_num_inc);


  for (unsigned int i = 0; i < _num_inc; ++i)
  {
    _n[i] = &adCoupledValue("incident", i);
    _E[i] = &adCoupledValue("incident_energies", i);
    _mu[i] = &getADMaterialProperty<Real>("mu" + (*getVar("incident", i)).name());
    _muE[i] = &getADMaterialProperty<Real>("mu" + (*getVar("incident_energies", i)).name());
    _sgn[i] = &getMaterialProperty<Real>("sgn" + (*getVar("incident", i)).name());
    _mass[i] = &getMaterialProperty<Real>("mass" + (*getVar("incident", i)).name());
  }
}

ADReal
HeatFluxKernel::computeQpResidual()
{
    // COMPUTE HEAT FLUX

  _Q_flux = 0.0;
  _v_th = 0.0;
  for (unsigned int i = 0; i < _num_inc; ++i)
  { 

    if ((*_sgn[i])[_qp] * -_grad_potential[_qp](0) > 0.0)
    {
      _a = 1.0;
    }
    else
    {
      _a = 0.0;
    }

    _v_th = std::sqrt(8 * 1.6E-19 * 2.0 / 3 * std::exp((*_E[i])[_qp] - (*_n[i])[_qp])
        / (M_PI * (*_mass[i])[_qp]));

    _Q_flux += _r_units * ((2. * _a - 1.) * (*_sgn[i])[_qp] * std::abs((*_muE[i])[_qp]) * -_grad_potential[_qp](0) 
    +5. / 6. * _v_th) * std::exp((*_E[i])[_qp]) * 6.022E23 * 1.6E-19;

  }
  return  -_test[_i][_qp]*_Q_flux;
}