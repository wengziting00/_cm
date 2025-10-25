import numpy as np

def solve_polynomial_roots(c: list):
    """
    使用數值方法求解多項式 f(x) = c[n]x^n + ... + c[1]x + c[0] 的根。

    參數:
    c (list): 多項式的係數列表，c[i] 是 x^i 的係數。
              注意: NumPy 要求係數是從最高次項開始排列。

    回傳:
    numpy.ndarray: 包含多項式根（可能為複數）的陣列。
    """
    
    # 檢查最高次項 c[n] 是否為 0，如果是，則需要移除尾部的 0
    while len(c) > 0 and c[-1] == 0:
        c.pop()
        
    if not c:
        # 如果所有係數都是 0 (即 f(x) = 0)
        print("警告：多項式 f(x) = 0，任何 x 都是解。")
        return np.array([])
        
    # 將係數列表反轉，以符合 NumPy 的要求：[最高次項, ..., 常數項]
    # 例如：x^3 + 2x + 1 的係數列表 c = [1, 2, 0, 1] 
    #       (c[0]=1, c[1]=2, c[2]=0, c[3]=1)
    #       反轉後變為: [1, 0, 2, 1]
    
    # 注意: 如果你的輸入 c 已經是 [c[n], ..., c[0]] 格式，則不需要反轉。
    # 根據你的描述「c[i] 代表 x^i 的係數」，所以我們需要反轉。
    reversed_c = c[::-1]
    
    print(f"原始係數 (c[i] for x^i): {c}")
    print(f"用於 NumPy 的係數 (從高次到低次): {reversed_c}")
    
    # 使用 numpy.roots 求解
    roots = np.roots(reversed_c)
    
    return roots
