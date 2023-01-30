//* This file is part of Zapdos, an open-source
//* application for the simulation of plasmas
//* https://github.com/shannon-lab/zapdos
//*
//* Zapdos is powered by the MOOSE Framework
//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADKernel.h"

class HeatFluxKernel : public ADKernel
{
public:

  static InputParameters validParams();
  
  HeatFluxKernel(const InputParameters & parameters);

protected:
  
  virtual ADReal computeQpResidual() override;
  
private:

  const Real _r_units;
  
  unsigned int _num_inc;
  std::vector<const ADVariableValue *> _n;
  std::vector<const ADVariableValue *> _E;
  std::vector<const ADMaterialProperty<Real> *> _mu;
  std::vector<const ADMaterialProperty<Real> *> _muE;
  std::vector<const MaterialProperty<Real> *> _sgn;
  std::vector<const MaterialProperty<Real> *> _mass;
  ADReal _a;
  ADReal _Q_flux;
  ADReal _v_th;
  const ADVariableGradient & _grad_potential;
};