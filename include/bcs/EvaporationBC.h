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

class EvaporationBC : public ADIntegratedBC
{
public: 
  static InputParameters validParams();
   
  EvaporationBC(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;
  
  const ADVariableValue & _T;
  const Real _r_units;
  const Real _a;
  const Real _b;
  const MaterialProperty<Real> & _m;
  const MaterialProperty<Real> & _kb;

  ADReal _evap_flux;

};