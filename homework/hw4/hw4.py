import numpy as np

def solve_polynomial_roots(c: list):

    while len(c) > 0 and c[-1] == 0:
        c.pop()
        
    if not c:
        # 如果所有係數都是 0 (即 f(x) = 0)
        print("警告：多項式 f(x) = 0，任何 x 都是解。")
        return np.array([])
        
    reversed_c = c[::-1]
    
    print(f"原始係數 (c[i] for x^i): {c}")
    print(f"用於 NumPy 的係數 (從高次到低次): {reversed_c}")
    
    roots = np.roots(reversed_c) 
    return roots

