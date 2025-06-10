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


def generate_bulk_free_energy_block(num_order_params_gr, num_order_params_regr, W=0.1, A=0.1, B=0.1, C=0.1, Sl=0.1):
    """
    生成MOOSE输入文件中的体相自由能块
    
    参数:
    num_order_params_gr: int, gr序参量的数量
    num_order_params_regr: int, regr序参量的数量
    W: float, 浓度项系数
    A: float, 序参量二次项系数
    B: float, 序参量四次项系数
    C: float, 序参量交叉项系数
    Sl: float, 耦合项系数
    """
    
    # 构建变量列表：c, gr0-grN, regr0-regrN
    variables = ['c']
    for i in range(num_order_params_gr):
        variables.append(f'gr{i}')
    for i in range(num_order_params_regr):
        variables.append(f'regr{i}')
    
    # 构建常数名称和表达式列表
    constant_names = ['W', 'A', 'B', 'C', 'Sl']
    constant_expressions = [str(W), str(A), str(B), str(C), str(Sl)]
    
    # 构建自由能表达式
    # 1. 浓度项
    f_conc = f"W*c^2*(1-c)^2"
    
    # 2. 序参量多项式项（所有gr和regr）
    f_order = []
    for i in range(num_order_params_gr):
        eta = f'gr{i}'
        f_order.append(f"(-A/2*{eta}^2+B/4*{eta}^4)")
    for i in range(num_order_params_regr):
        eta = f'regr{i}'
        f_order.append(f"(-A/2*{eta}^2+B/4*{eta}^4)")
    f_order_sum = "+".join(f_order)
    
    # 3. 序参量交叉项（所有序参量之间的交叉项）
    f_cross = []
    # gr之间的交叉
    for i in range(num_order_params_gr):
        for j in range(num_order_params_gr):
            if i != j:
                f_cross.append(f"gr{i}^2*gr{j}^2")
    # regr之间的交叉
    for i in range(num_order_params_regr):
        for j in range(num_order_params_regr):
            if i != j:
                f_cross.append(f"regr{i}^2*regr{j}^2")
    # gr与regr之间的交叉
    for i in range(num_order_params_gr):
        for j in range(num_order_params_regr):
            f_cross.append(f"gr{i}^2*regr{j}^2")
            f_cross.append(f"regr{j}^2*gr{i}^2")
    
    f_cross_sum = f"C*("+"+".join(f_cross)+")" if f_cross else "0"
    
    # 4. 耦合项（所有序参量的平方项）
    eta_terms = []
    for i in range(num_order_params_gr):
        eta_terms.append(f"gr{i}^2")
    for i in range(num_order_params_regr):
        eta_terms.append(f"regr{i}^2")
    eta_sum = "+".join(eta_terms) if eta_terms else "0"
    f_coupling = f"1/4+Sl*c*(1-({eta_sum}))"
    
    # 组合完整表达式
    full_expression = f"{f_conc}+{f_order_sum}+{f_cross_sum}+{f_coupling}"
    
    # 生成MOOSE块
    moose_block = f"""
  [./bulk_free_energy]
    type = DerivativeParsedMaterial
    property_name = f_bulk
    coupled_variables = '{" ".join(variables)}'
    constant_names = '{" ".join(constant_names)}'
    constant_expressions = '{" ".join(constant_expressions)}'
    expression = {full_expression}
    derivative_order = 2
    outputs = exodus
  [../]
"""
    
    return moose_block

def generate_f_local(num_order_params_gr, num_order_params_regr):
    variables = ['c']
    variablesbulk = ['c']
    variablesstored = []
    for i in range(num_order_params_gr):
        variables.append(f'gr{i}')
        variablesbulk.append(f',gr{i}')
        variablesstored.append(f'gr{i},')
    for i in range(num_order_params_regr):
        variables.append(f'regr{i}')
        variablesbulk.append(f',regr{i}')
        if i == num_order_params_regr - 1:
            variablesstored.append(f'regr{i}')
        else:
            variablesstored.append(f'regr{i},')

    #去掉variablesstored中的最后一个逗号,注意regr{i},是一体的


    moose_block = f"""
  [./f_local]
    type = DerivativeParsedMaterial
    property_name = f_local
    coupled_variables = '{" ".join(variables)}'
    expression = 'f_bulk + f_stored'
    material_property_names = 'f_bulk({"".join(variablesbulk)}) f_stored({"".join(variablesstored)})'
    outputs = exodus
  [../]
"""
    return moose_block

# 使用示例
if __name__ == "__main__":
    # 生成每个序参量有独立位错密度的版本
    num_order_params_gr = 3
    num_order_params_regr = 1
    # 生成共享位错密度的版本
    print("=== 共享位错密度版本 ===")
    moose_block1 = generate_f_stored_block_shared_rho(
        num_order_params_gr=num_order_params_gr,
        num_order_params_regr=num_order_params_regr,
        G='${shear_modulus}',
        b='${burgers_vector}'
    )
    print(moose_block1)
    moose_block2 = generate_bulk_free_energy_block(
    num_order_params_gr=num_order_params_gr,
    num_order_params_regr=num_order_params_regr,
    W=0.1,
    A=0.2,
    B=0.3,
    C=0.4,
    Sl=0.5
    )
    print(moose_block2)
    moose_block3 = generate_f_local(
        num_order_params_gr=num_order_params_gr,
        num_order_params_regr=num_order_params_regr
    )
    print(moose_block3)