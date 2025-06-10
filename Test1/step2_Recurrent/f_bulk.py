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

# 使用示例
if __name__ == "__main__":
    # 可以分别指定gr和regr的数量
    moose_block = generate_bulk_free_energy_block(
        num_order_params_gr=3,
        num_order_params_regr=2,
        W=0.1,
        A=0.2,
        B=0.3,
        C=0.4,
        Sl=0.5
    )
    print(moose_block)
