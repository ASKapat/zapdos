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

#include "ADIntegratedBC.h"

class HagelaarEnergyBC_HS : public ADIntegratedBC
{
public:
  static InputParameters validParams();

  HagelaarEnergyBC_HS(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const Real _r_units;
  const Real & _r;

  const ADVariableGradient & _grad_potential;
  const ADVariableValue & _n_HS;
  const ADMaterialProperty<Real> & _muip;

  const MaterialProperty<Real> & _massip;
  const MaterialProperty<Real> & _sgnip;
  const MaterialProperty<Real> & _e;

  const ADMaterialProperty<Real> & _muE_HS;

  Real _a;
  ADReal _v_thermal;
};
