#include "MyADBulkEnergyD2.h"

registerMooseObject("MissCatApp", MyADBulkEnergyD2);

InputParameters
MyADBulkEnergyD2::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Computes bulk free energy density with automatic differentiation for an arbitrary number of grain order parameters.");

  // Require a vector of coupled variables for grain order parameters
  params.addRequiredCoupledVar("grs", "Vector of grain order parameters (e.g., 'gr0 gr1 gr2 ...')");

  // Material properties (can be coupled from other materials or given as constants)
  params.addParam<Real>("mu", "Material property mu");
  params.addParam<Real>("gamma_asymm", "Material property gamma_asymm");

  params.addParam<MaterialPropertyName>("f_name", "f_bulk", "Material property name for the bulk free energy density");

  return params;
}

MyADBulkEnergyD2::MyADBulkEnergyD2(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _num_grs(coupledComponents("grs")),
    _gr_vars(_num_grs),
    _gr_var_names(_num_grs),
    _mu(getParam<Real>("mu")),
    _gamma_asymm(getParam<Real>("gamma_asymm")),
    // 使用 declareADProperty 声明 AD 材料属性
    _f_bulk(declareADProperty<Real>(getParam<MaterialPropertyName>("f_name"))),
    // 声明导数属性
    _df_bulk_dgr(_num_grs)
{
  // 存储耦合变量和变量名
  for (unsigned int i = 0; i < _num_grs; ++i)
  {
    _gr_vars[i] = &adCoupledValue("grs", i); // 使用 adCoupledValue 获取 AD 变量值
    _gr_var_names[i] = getVar("grs", i)->name();
    
    // 为每个序参数声明导数属性
    _df_bulk_dgr[i] = &declareADProperty<Real>(derivativePropertyNameFirst(getParam<MaterialPropertyName>("f_name"), _gr_var_names[i]));
  }
  // 不需要手动声明导数属性，AD 框架会自动处理
}

void
MyADBulkEnergyD2::computeQpProperties()
{
  ADReal sum_gr_sq = 0.0;
  ADReal sum_gr_pow4 = 0.0;
  for (unsigned int i = 0; i < _num_grs; ++i)
  {
    ADReal gr_i = (*_gr_vars[i])[_qp];
    sum_gr_sq += gr_i * gr_i;
    sum_gr_pow4 += std::pow(gr_i, 4);
  }

  ADReal sum_gr_pair_interaction = 0.0;
  for (unsigned int i = 0; i < _num_grs; ++i)
  {
    for (unsigned int j = i + 1; j < _num_grs; ++j) // j > i to avoid double counting and self-interaction in this term
    {
      ADReal gr_i_sq = std::pow((*_gr_vars[i])[_qp], 2);
      ADReal gr_j_sq = std::pow((*_gr_vars[j])[_qp], 2);
      sum_gr_pair_interaction += gr_i_sq * gr_j_sq;
    }
  }

  _f_bulk[_qp] = _mu * (-0.5 * sum_gr_sq + 0.25 * sum_gr_pow4 + _gamma_asymm * sum_gr_pair_interaction) + 0.25;
  
  // 计算导数属性 - 在AD框架中，这些导数会自动计算
  // 我们不需要手动计算导数，但需要将导数值赋给声明的导数属性
  for (unsigned int i = 0; i < _num_grs; ++i)
  {
    // 在AD框架中，导数会自动计算并传播，我们只需要将_f_bulk的导数赋值给对应的导数属性
    // 这里不需要手动计算导数，AD框架会自动处理
    // 但我们需要确保导数属性被正确声明和赋值
    (*_df_bulk_dgr[i])[_qp] = _f_bulk[_qp].derivatives()[i];
  }
}
