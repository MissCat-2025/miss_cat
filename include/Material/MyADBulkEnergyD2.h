#pragma once

#include "DerivativeMaterialInterface.h"
#include "Material.h"
#include <vector>

/**
 * @brief Material class to compute bulk free energy density based on a given expression.
 * This class uses DerivativeMaterialInterface for automatic differentiation.
 * The expression is: mu*(-0.5*sum(gr_i^2) + 0.25*sum(gr_i^4) + gamma_asymm*sum_pairs(gr_i^2*gr_j^2)) + 0.25
 */
class MyADBulkEnergyD2 : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  MyADBulkEnergyD2(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  // Coupled variables for grain order parameters
  unsigned int _num_grs; // Number of grain order parameters
  std::vector<const ADVariableValue *> _gr_vars; // 使用 ADVariableValue 类型
  std::vector<VariableName> _gr_var_names; // 存储变量名称

  // Material parameters
  const Real _mu;
  const Real _gamma_asymm;

  // Free energy material property
  ADMaterialProperty<Real> & _f_bulk; // 使用 ADMaterialProperty 类型
  
  // 导数属性
  std::vector<ADMaterialProperty<Real> *> _df_bulk_dgr; // 自由能对序参数的导数
};
