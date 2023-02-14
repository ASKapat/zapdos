#include "SputterBC.h"

registerADMooseObject("ZapdosApp", SputterBC);

InputParameters
SputterBC::validParams()
{
InputParameters params = ADIntegratedBC::validParams();
  params.addRequiredCoupledVar("incident", "Ion species(es) incident on surface");
  params.addRequiredCoupledVar("incident_energy", "Energy of incident ion(s)");
  params.addRequiredCoupledVar("sputtered", "Surface species sputtered");
  params.addRequiredParam<Real>("sputter_Z", "Atomic number of sputtered species");
  params.addRequiredCoupledVar("potential", "The electric potential");
  params.addRequiredParam<Real>("position_units", "Units of position.");
  params.addRequiredParam<std::vector<Real>>("incident_Z", "Atomic number(s) of incident species");
  params.addRequiredParam<std::vector<Real>>("r_ion", "reflection coefficient(s) of ions");
  params.addRequiredParam<std::vector<Real>>("q", "Free parameter for Eckstein sputtering.");
  params.addRequiredParam<std::vector<Real>>("Threshold", "Threshold parameter for Eckstein sputtering");
  params.addRequiredParam<std::vector<Real>>("mu", "Eckstein sputtering empirical exponent");
  params.addRequiredParam<std::vector<Real>>("lambda", "Eckstein sputtering third 4th parameter");
  params.addRequiredParam<std::vector<Real>>("Sputter_fraction","fraction of total sputtering yield");
  
  return params;
}

SputterBC::SputterBC(const InputParameters & parameters)
	: ADIntegratedBC(parameters),
	  _r_units(1. /getParam<Real>("position_units")),
	  _kb(getMaterialProperty<Real>("k_boltz")),
	  _Zs(getParam<Real>("sputter_Z")),
	  
	  
	 

	  
	  
	  // Coupled Variables
	  _grad_potential(adCoupledGradient("potential")),
	  _e(getMaterialProperty<Real>("e")),
	  _N_A(getMaterialProperty<Real>("N_A")),
	  _n_s(adCoupledValue("sputtered")),
	  _m_s(getMaterialProperty<Real>("mass"+(*getVar("sputtered", 0)).name())),
	  _sgn_s(getMaterialProperty<Real>("sgn"+(*getVar("sputtered", 0)).name())),
	  
	  
	  _r_ion(getParam<std::vector<Real>>("r_ion")),
	  _Zi(getParam<std::vector<Real>>("incident_Z")),
	  _q(getParam<std::vector<Real>>("q")),
	  _Eth(getParam<std::vector<Real>>("Threshold")),
	  _mu2(getParam<std::vector<Real>>("mu")),
	  _lamb2(getParam<std::vector<Real>>("lambda")),
	  _sput_frac(getParam<std::vector<Real>>("Sputter_fraction"))
	  
	
	  
	  /*Need to get sputtering coefficient from get material property
	  _se_coeff(getMaterialProperty<Real>("se_coeff"))
	  */
{
  _num_ions = coupledComponents("incident");
  _n_ip.resize(_num_ions);
  _Eip.resize(_num_ions);
  _muip.resize(_num_ions);
  _massip.resize(_num_ions);
  _sgnip.resize(_num_ions);
//   _E_TF.resize(_num_ions);

  
//   std::cout << _num_ions << std::endl;
  for (unsigned int i = 0; i < _num_ions; ++i)
  {
    _n_ip[i] = &adCoupledValue("incident", i);
    _Eip[i] = &adCoupledValue("incident_energy", i);
    _muip[i] = &getADMaterialProperty<Real>("mu" + (*getVar("incident", i)).name());
    _massip[i] = &getMaterialProperty<Real>("mass" + (*getVar("incident", i)).name());
    _sgnip[i] = &getMaterialProperty<Real>("sgn" + (*getVar("incident", i)).name());
//     _E_TF[i] = 30.74*(_m_s[_qp]+(*_massip[i])[_qp])/(_m_s[_qp])*_Zs*_Zi[i]*std::sqrt(std::pow(_Zs, 2/3)+std::pow(_Zi[i], 2/3));
//      std::cout << i << std::endl;

  }
  
  
}

ADReal
SputterBC::computeQpResidual()
{
  if (_normals[_qp] * 1.0 * -_grad_potential[_qp] > 0.0)
    {
      _a = 1.0;
    }
    else
    {
      _a = 0.0;
    }
  _sputt_flux = 0;
  _ion_flux = 0;
  _Y = 0;
//   std::cout << _Y<< std::endl;

  for (unsigned int i = 0; i < _num_ions; ++i)
  { 
//     std::cout << i << std::endl;

  _E_TF = 30.74*((_m_s[_qp]+(*_massip[i])[_qp])/(_m_s[_qp]))*_Zs*_Zi[i]*std::sqrt(std::pow(_Zs, 2/3)+std::pow(_Zi[i], 2/3));
  
  _v_th_i = std::sqrt(8 * _e[_qp] * 2.0 / 3 * std::exp((*_Eip[i])[_qp] - (*_n_ip[i])[_qp])
        / (M_PI * (*_massip[i])[_qp]));

  _ion_flux = (1-_r_ion[i])/(1+_r_ion[i])* std::exp((*_n_ip[i])[_qp]) * (0.5 * _v_th_i + _a  *
         ((*_sgnip[i])[_qp] * (*_muip[i])[_qp] * -_grad_potential[_qp] *_r_units * _normals[_qp]));
         

  _red_e = (std::exp((*_Eip[i])[_qp] - (*_n_ip[i])[_qp]))/(_E_TF);
  
  _s_n = 0.5*std::log(1+1.2288*_red_e)/(_red_e + 0.1728*std::sqrt(_red_e)+0.008*std::pow(_red_e, 0.1504));

//   _Y = (_q[i]*_s_n*std::pow(std::exp((*_Eip[i])[_qp] - (*_n_ip[i])[_qp])/(_Eth[i])-1), (_mu2[i]))/
//    ((_lamb2[i])+std::pow(std::exp((*_Eip[i])[_qp] - (*_n_ip[i])[_qp])/(_Eth[i])-1, (_mu2[i])));
   
   _Y = (_q[i]*_s_n*std::pow(std::exp((*_Eip[i])[_qp] - (*_n_ip[i])[_qp])/(_Eth[i])-1,_mu2[i]))/
   (_lamb2[i]+std::pow(std::exp((*_Eip[i])[_qp] - (*_n_ip[i])[_qp])/(_Eth[i])-1,_mu2[i]));
   
   std::cout << _Y<< std::endl;
   std::cout << _red_e<< std::endl;
   std::cout << _s_n<< std::endl;
   std::cout << std::exp((*_Eip[i])[_qp] - (*_n_ip[i])[_qp])<<std::endl;
  _sputt_flux += _ion_flux*_Y*_sput_frac[i];
  }
//   std::cout << _r_units * (1. - _a * _sgn_s[_qp]) * _sputt_flux << std::endl;
//   std::cout <<  _ion_flux << std::endl;
//   std::cout << _Y<< std::endl;
//   std::cout << _sputt_flux << std::endl;
  return -_test[_i][_qp] * _r_units * (1. - _a * _sgn_s[_qp]) * _sputt_flux;
  

  
}