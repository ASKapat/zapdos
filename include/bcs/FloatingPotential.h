#pragma once

#include "ADIntegratedBC.h"

/**
 * Boundary condition of a Dirichlet type
 *
 * Sets the value at the node to the value of a Postprocessor
 */
class FloatingPotential : public ADIntegratedBC
{
public:
  static InputParameters validParams();

  FloatingPotential(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  /// The value for this BC
  
  const Real _Vp;
  const VariableValue & _em;
  const VariableValue & _Te;
  
  const std::vector<Real> & _massip;
  std::string _potential_units;
  const Real _r_units;
  unsigned int _num_ions;
  std::vector<const VariableValue *> _n_ip;
  std::vector<const VariableValue *> _Eip;
  Real _voltage_scaling;
  Real _p;
  Real _ion_flux;
  Real _ion_flux_tot;
  Real _V_f;
  
};