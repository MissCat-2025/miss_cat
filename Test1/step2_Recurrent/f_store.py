

def generate_f_stored_block_shared_rho(num_order_params_gr, num_order_params_regr, G='${shear_modulus}', b='${burgers_vector}'):
    """
    生成MOOSE输入文件中的储存能材料块（共享位错密度版本）
    
    参数:
    num_order_params_gr: int, gr序参量的数量
    num_order_params_regr: int, regr序参量的数量
    G: str, 剪切模量（可以是变量或数值）
    b: str, 柏格斯矢量（可以是变量或数值）
    """
    
    # 构建变量列表：gr0-grN, regr0-regrN
    variables = []
    for i in range(num_order_params_gr):
        variables.append(f'gr{i}')
    for i in range(num_order_params_regr):
        variables.append(f'regr{i}')
    
    # 构建常数名称和表达式列表
    constant_names = ['G', 'b_v2']
    constant_expressions = [str(G), str(b)]
    
    # 构建储存能表达式（所有序参量共享一个rho）
    eta_squared_terms = []
    
    # gr项
    for i in range(num_order_params_gr):
        eta_squared_terms.append(f"gr{i}^2")
    
    # regr项
    for i in range(num_order_params_regr):
        eta_squared_terms.append(f"regr{i}^2")
    
    # 组合完整表达式
    eta_sum = "+".join(eta_squared_terms) if eta_squared_terms else "0"
    full_expression = f"0.5*G*b_v2*rho*({eta_sum})"
    
    # 生成MOOSE块
    moose_block = f"""
  [./f_stored]
    type = DerivativeParsedMaterial
    property_name = f_stored
    coupled_variables = '{" ".join(variables)}'
    material_property_names = 'rho'
    constant_names = '{" ".join(constant_names)}'
    constant_expressions = '{" ".join(constant_expressions)}'
    expression = {full_expression}
    derivative_order = 2
    outputs = exodus
  [../]
"""
    
    return moose_block

# 使用示例
if __name__ == "__main__":
    # 生成每个序参量有独立位错密度的版本
    
    # 生成共享位错密度的版本
    print("=== 共享位错密度版本 ===")
    moose_block2 = generate_f_stored_block_shared_rho(
        num_order_params_gr=3,
        num_order_params_regr=2,
        G='${shear_modulus}',
        b='${burgers_vector}'
    )
    print(moose_block2)
