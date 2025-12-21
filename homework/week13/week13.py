import numpy as np
from collections import Counter

def solve_ode_general(coefficients):
    roots = np.roots(coefficients)
    rounded_roots = []
    for r in roots:
        real_part = round(r.real, 10)
        imag_part = round(r.imag, 10)
        if abs(imag_part) < 1e-10:
            rounded_roots.append(real_part + 0j)
        else:
            rounded_roots.append(complex(real_part, imag_part))
    
    root_counts = Counter(rounded_roots)
    
    solutions = []
    c_idx = 1
    processed_complex = set()

    unique_roots = sorted(root_counts.keys(), key=lambda x: (x.real, x.imag))
    for r in unique_roots:
        if r in processed_complex:
            continue
        
        m = root_counts[r] 

        if abs(r.imag) > 1e-10:
            alpha = r.real
            beta = abs(r.imag)
            conjugate_root = complex(r.real, -r.imag)
            processed_complex.add(conjugate_root)
            
            for k in range(m):
                x_term = f"x^{k}" if k > 0 else ""
                e_term = f"e^({alpha}x)" if alpha != 0 else ""
                solutions.append(f"C_{c_idx}{x_term}{e_term}cos({beta}x)")
                c_idx += 1
                solutions.append(f"C_{c_idx}{x_term}{e_term}sin({beta}x)")
                c_idx += 1
        else:
            val = r.real
            for k in range(m):
                x_term = f"x^{k}" if k > 0 else ""
                if val == 0:
                    solutions.append(f"C_{c_idx}{x_term}")
                else:
                    solutions.append(f"C_{c_idx}{x_term}e^({val}x)")
                c_idx += 1
    final_str = " + ".join(solutions)
    final_str = final_str.replace("x^1e", "xe").replace("x^1cos", "xcos").replace("x^1sin", "xsin")
    return f"y(x) = {final_str}"

# --- 以下為測試主程式 ---

# 範例測試 (1): 實數單根: y'' - 3y' + 2y = 0
print("--- 實數單根範例 ---")
coeffs1 = [1, -3, 2]
print(f"方程係數: {coeffs1}")
print(solve_ode_general(coeffs1))

# 範例測試 (2): 實數重根: y'' - 4y' + 4y = 0
print("\n--- 實數重根範例 ---")
coeffs2 = [1, -4, 4]
print(f"方程係數: {coeffs2}")
print(solve_ode_general(coeffs2))

# 範例測試 (3): 複數共軛根: y'' + 4y = 0
print("\n--- 複數共軛根範例 ---")
coeffs3 = [1, 0, 4]
print(f"方程係數: {coeffs3}")
print(solve_ode_general(coeffs3))

# 範例測試 (4): 複數重根 (二重): (D^2 + 1)^2 y = 0
print("\n--- 複數重根範例 ---")
coeffs4 = [1, 0, 2, 0, 1]
print(f"方程係數: {coeffs4}")
print(solve_ode_general(coeffs4))

# 範例測試 (5): 高階重根: y''' - 6y'' + 12y' - 8y = 0
print("\n--- 高階重根範例 ---")
coeffs5 = [1, -6, 12, -8]
print(f"方程係數: {coeffs5}")
print(solve_ode_general(coeffs5))
