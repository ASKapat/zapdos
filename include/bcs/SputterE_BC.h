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

class SputterE_BC : public ADIntegratedBC
{
public:
  static InputParameters validParams();

  SputterE_BC(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const Real _r_units;
  
  const MaterialProperty<Real> & _kb;
  const Real & _Zs;
  const Real & _E_sb;
  const ADVariableGradient & _grad_potential;
  const MaterialProperty<Real> & _e;
  const MaterialProperty<Real> & _N_A;
  const ADVariableValue & _n_s;
  const MaterialProperty<Real> & _m_s;
  const MaterialProperty<Real> & _sgn_s;
  unsigned int _num_ions;
  
  const std::vector<Real> & _r_ion;
  const std::vector<Real> & _Zi;
  const std::vector<Real> & _q;
  const std::vector<Real> & _Eth;
  const std::vector<Real> & _mu2;
  const std::vector<Real> & _lamb2;
  const std::vector<Real> & _sput_frac;
  
  std::vector<const ADVariableValue *> _n_ip;
  std::vector<const ADVariableValue *> _Eip;
  std::vector<const ADMaterialProperty<Real> *> _muip;
  std::vector<const MaterialProperty<Real> *> _massip;
  std::vector<const MaterialProperty<Real> *> _sgnip;
//   std::vector<ADReal> _E_TF;
  
  
  
  Real _a;
  ADReal _sputt_flux;
  ADReal _ion_flux;
  
  ADReal _E_TF;
  ADReal _v_th_i;
  
  ADReal _red_e;
  ADReal _s_n;
  ADReal _Y;
  
  };