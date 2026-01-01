import numpy as np

def solve_linear_equations():
    # 定義係數矩陣 A
    A = np.array([
        [2, 1, -1],
        [-3, -1, 2],
        [-2, 1, 2]
    ])

    # 定義常數向量 b
    b = np.array([8, -11, -3])

    try:
        # 使用 numpy 的線性代數求解器
        x = np.linalg.solve(A, b)
        
        print("求解成功！結果如下：")
        print(f"x = {x[0]:.2f}")
        print(f"y = {x[1]:.2f}")
        print(f"z = {x[2]:.2f}")
        
        # 驗證結果：計算 Ax 是否等於 b
        verification = np.allclose(np.dot(A, x), b)
        print(f"\n驗證 (Ax = b): {'通過' if verification else '失敗'}")
        
    except np.linalg.LinAlgError:
        print("錯誤：矩陣為奇異矩陣 (Singular Matrix)，方程組有無窮多解或無解。")

if __name__ == "__main__":
    solve_linear_equations()
