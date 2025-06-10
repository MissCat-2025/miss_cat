import numpy as np
import matplotlib.pyplot as plt

# 设置中文显示
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']
plt.rcParams['axes.unicode_minus'] = False

# 物理常数
k_B = 8.617333262e-5  # 玻尔兹曼常数 (eV/K)

# 已知参数 (全部转换为数值)
L = 3.0e-8       # Distance between planes (m)
a = 3.42e-10     # Lattice constant (m)
B_i = 5.0e-3     # Relative bias between interstitials and vacancies
D_i0 = 3.81e-9   # Interstitial diffusivity pre-factor (m²/s)
D_v0 = 5.72e-8   # Vacancy diffusivity pre-factor (m²/s)
epsilon_i = 0.27 # Interstitial migration enthalpy (eV)
epsilon_v = 1.08 # Vacancy migration enthalpy (eV)
r_iv = 8.25e-9   # Radius of recombination volume (m)
Omega = a**3/2   # Atomic volume (m³)



def calculate_dislocation_density(T, J_dot):
    """计算位错密度
    T: 温度 (K)
    J_dot: 裂变率 (fissions/cm³/s)
    返回: 位错密度 (m⁻²)
    ρ_d=(3/La)^(2/3) (α_i/(B_i D_i ))^(1/3) ((KD_v)/α_r )^(1/6),
    """
    # 计算K值
    
    # 计算各个常数项
    term1 = (3/(L*a))**(2/3)
    term2 = ((np.sqrt(2)/a**2)/(B_i))**(1/3)
    #温度项
    term3 = D_v0**(1/6)*D_i0**(-1/3)
    KKT = 1/(6*k_B)*(2*epsilon_i-epsilon_v)
    term4 = np.exp(KKT/T)
    print("\K系数为:")
    print(KKT)

    term5 = (1 / (5.0e23))**(1/6)
    term6 = (J_dot)**(1/6)
    term7 = (Omega/(4*np.pi*r_iv))**(1/6)
    rho_d1 = term1 * term2 * term3 * term5*term7
    print("\常数系数为:")
    print(rho_d1)
    # 计算最终结果
    rho_d = term1 * term2 * term3 * term4 * term5 * term6*term7
    
    # 打印各项的值
    print(f"最终位错密度 ρ_d = {rho_d:.2e} m⁻²")
    #变为纳米单位：
    rho_d2 = rho_d*1e-18
    print(f"最终位错密度 ρ_d = {rho_d2:.2e} nm⁻²")
    
    return rho_d

T = 373
J_dot = 0.98e13
calculate_dislocation_density(T, J_dot)
