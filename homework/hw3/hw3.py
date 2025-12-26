import cmath

def root3(a, b, c, d):
    # 將三次方程式 ax^3 + bx^2 + cx + d = 0 轉換為簡化形式 t^3 + pt + q = 0
    if a == 0:
        raise ValueError("a 不能為 0")

    # 除以 a，使首項係數為 1
    b /= a
    c /= a
    d /= a

    # 簡化
    p = c - b**2 / 3
    q = (2 * b**3) / 27 - (b * c) / 3 + d

    # 計算判別式 Δ
    # Δ > 0：一個實根
    # Δ = 0：三個實根
    # Δ < 0: 三個不同的實根
    Δ = (q / 2)**2 + (p / 3)**3

    # 使用 Cardano 公式
    u_cube = -q / 2 + cmath.sqrt(Δ)
    v_cube = -q / 2 - cmath.sqrt(Δ)

    u = u_cube ** (1/3) if u_cube != 0 else 0
    v = v_cube ** (1/3) if v_cube != 0 else 0

    def cube_root(z):
        return z ** (1/3) if z.imag == 0 else cmath.exp(cmath.log(z) / 3)

    u = cube_root(u_cube)
    v = cube_root(v_cube)

    # 三個根（使用立方根單位根）
    omega = complex(-0.5, cmath.sqrt(3)/2)  # ω = e^(2πi/3)

    t1 = u + v
    t2 = u * omega + v * omega**2
    t3 = u * omega**2 + v * omega

    # 回代 x = t - b/3
    x1 = t1 - b / 3
    x2 = t2 - b / 3
    x3 = t3 - b / 3

    return (x1, x2, x3)
