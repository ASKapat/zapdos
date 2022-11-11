//* This file is part of Zapdos, an open-source
//* application for the simulation of plasmas
//* https://urldefense.com/v3/__https://github.com/shannon-lab/zapdos__;!!DZ3fjg!sH2ulMsTKshznG7u93hGEbXsjzERnqjXaL6VLUIMdNAAbzTmOnf4rRx3CG56APwVZA$ 
//*
//* Zapdos is powered by the MOOSE Framework
//* https://urldefense.com/v3/__https://www.mooseframework.org__;!!DZ3fjg!sH2ulMsTKshznG7u93hGEbXsjzERnqjXaL6VLUIMdNAAbzTmOnf4rRx3CG6mxSP5BA$ 
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://urldefense.com/v3/__https://www.gnu.org/licenses/lgpl-2.1.html__;!!DZ3fjg!sH2ulMsTKshznG7u93hGEbXsjzERnqjXaL6VLUIMdNAAbzTmOnf4rRx3CG6Yp_I17g$ 

#pragma once

#include "Material.h"
/* #include "LinearInterpolation.h" */
#include "SplineInterpolation.h"

template <bool is_ad>
class HeavySpecies_ETempl : public Material
{
public:
  HeavySpecies_ETempl(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual void computeQpProperties() override;

  Real _user_massHeavy;
  Real _user_sgnHeavy;

  std::string _potential_units;
  MaterialProperty<Real> & _massHeavy;
  GenericMaterialProperty<Real, is_ad> & _temperatureHeavy;
  MaterialProperty<Real> & _sgnHeavy; 
    
  GenericMaterialProperty<Real, is_ad> & _muHeavy;          // Replaces _muArp
  GenericMaterialProperty<Real, is_ad> & _diffHeavy;        // Replaces _diffArp
  GenericMaterialProperty<Real, is_ad> & _diffHeavy_E;      // HS  energy diff coeff
  GenericMaterialProperty<Real, is_ad> & _muHeavy_E;        // HS energy mu coeff
  MaterialProperty<Real> & _sgnHeavy_E;                
  Real _voltage_scaling;
                    
  

  const VariableValue & _E_heavy;
  const VariableValue & _heavy_n;




  Real _time_units;
  bool _calc_mobility;
  bool _calc_diffusivity;


};

typedef HeavySpecies_ETempl<false> HeavySpecies_E;
typedef HeavySpecies_ETempl<true> ADHeavySpecies_E;