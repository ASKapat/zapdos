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

class FloatingPotentialNeumann : public ADIntegratedBC
{
public:
  static InputParameters validParams();

  FloatingPotentialNeumann(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const ADVariableValue & _em;
  const ADVariableValue & _mean_en;
  const std::string & _potential_units;
  const Real _r_units;
  const MaterialProperty<Real> & _eps;
  const MaterialProperty<Real> & _N_A;
  const ADMaterialProperty<Real> & _muem;
  const MaterialProperty<Real> & _e;
  const MaterialProperty<Real> & _massem;
  const MaterialProperty<Real> & _kb;
  ADReal _a_e;
  Real _voltage_scaling;
  Real _num_ions;
  std::vector<const ADVariableValue *> _ip;
  std::vector<const ADVariableValue *> _Eip;
  std::vector<const ADMaterialProperty<Real> *> _muip;
  std::vector<const MaterialProperty<Real> *> _massip;
  std::vector<const MaterialProperty<Real> *> _sgnip;

  
  std::vector<ADReal> _a;
  ADReal _ion_therm_drift;
  ADReal _advect;
  ADReal _v_th_i;
  ADReal _v_e_th;
  
};