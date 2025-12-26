import cmath

def root2(a, b, c):
    discriminant = cmath.sqrt(b**2 - 4*a*c)# 計算平方根的部分

    root1 = (-b + discriminant) / (2*a)# 第一
    root2 = (-b - discriminant) / (2*a)# 第二

    # 把 root1 代入 f(x)
    f1 = a*root1**2 + b*root1 + c
    f2 = a*root2**2 + b*root2 + c

    # 用 cmath.isclose 判斷是不是接近 0
    if not cmath.isclose(f1, 0, rel_tol=1e-9, abs_tol=0.0):
        print("Warning: root1 代入後不接近 0")
    if not cmath.isclose(f2, 0, rel_tol=1e-9, abs_tol=0.0):
        print("Warning: root2 代入後不接近 0")

    return root1, root2
