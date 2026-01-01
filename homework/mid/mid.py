import numpy as np
from scipy.optimize import fsolve

def equations(vars):
    x, y = vars
    # 將方程式寫成「等於 0」的形式
    eq1 = x**2 + y - 10
    eq2 = x + y**2 - 7
    return [eq1, eq2]

def solve_nonlinear_equations():
    # 非線性求解需要一個「初始猜測值」 (initial guess)
    initial_guess = [1, 1]
    
    try:
        # 使用 fsolve 求解
        result = fsolve(equations, initial_guess)
        
        print("求解成功 結果如下：")
        print(f"x = {result[0]:.4f}")
        print(f"y = {result[1]:.4f}")
        
        # 驗證代入原方程是否接近於 0
        test_res = equations(result)
        print(f"\n驗證 (是否接近0): {np.allclose(test_res, [0, 0])}")
        
    except Exception as e:
        print(f"求解失敗: {e}")

if __name__ == "__main__":
    solve_nonlinear_equations()
