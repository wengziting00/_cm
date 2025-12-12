# 離散傅立葉變換 (DFT) 與逆變換 (IDFT) 實現及可逆性驗證

這個專案使用 Python 和 NumPy 庫，從頭實現了離散傅立葉變換 (DFT) 及其逆變換 (IDFT)。核心目的是透過矩陣運算來高效計算這兩種變換，並驗證它們互為逆運算 (即 $IDFT(DFT(x)) \approx x$) 的數學性質。

1. 數學原理

A. 離散傅立葉變換 (DFT)

DFT 將一個長度為 $N$ 的時域離散信號 $x[n]$ 轉換為頻域信號 $X[k]$。

$$X[k] = \sum_{n=0}^{N-1} x[n] \cdot e^{-i 2 \pi k n / N}$$

B. 離散傅立葉逆變換 (IDFT)

IDFT 將頻域信號 $X[k]$ 逆轉換回原始時域信號 $x[n]$。

$$x[n] = \frac{1}{N} \sum_{k=0}^{N-1} X[k] \cdot e^{+i 2 \pi k n / N}$$

 2. 程式碼實現 (使用 NumPy 矩陣運算)

本實現避免了多重迴圈，利用 NumPy 的向量化操作和矩陣乘法 (`.dot()`) 來計算變換矩陣 $M$ 與信號向量的乘積。

#

2.1. DFT 函數 (`dft`)

在 DFT 中，我們定義了變換矩陣 $M$ 的元素 $M_{k,n} = e^{-i 2 \pi k n / N}$。

```python
import numpy as np

def dft(x):
    N = len(x)
    n = np.arange(N)
    k = n.reshape((N, 1)) # 使 k 成為 (N, 1) 的列向量

    # 核心變換矩陣 M: M_kn = exp(-2j * pi * k * n / N)
    M = np.exp(-2j * np.pi * k * n / N)

    # 矩陣乘法: X = M * x (N x N 矩陣 . N x 1 向量)
    X = M.dot(x)

    return X


2.2. IDFT 函數 (idft)在 IDFT 中，變換矩陣的指數符號相反 (變為正號)，並且結果需要除以信號長度 $N$。

def idft(X): #
    N = len(X)
    k = np.arange(N)
    n = k.reshape((N, 1))

    # 核心變換矩陣 M': M'_nk = exp(+2j * pi * k * n / N)
    M = np.exp(2j * np.pi * k * n / N)

    # 矩陣乘法: x_unscaled = M' * X
    x_unscaled = M.dot(X)

    # 縮放因子 (1/N)
    x = x_unscaled / N

    return x


 3. 可逆性驗證verify_dft_inverse 函數執行 $DFT \to IDFT$ 的過程，並計算原始信號與恢復信號之間的 L2 誤差。

3.1. 驗證函數 (verify_dft_inverse)

def verify_dft_inverse(f):
    print("\n--- 驗證 DFT/IDFT 可逆性 ---")
    
    # 1. DFT
    F = dft(f)
    print(f"原始信號長度 N={len(f)}")
    
    # 2. IDFT
    f_recovered = idft(F)
    f_recovered_real = f_recovered.real # 取實部 (原始信號是實數)
    
    # 3. 計算 L2 誤差 (歐幾里得距離)
    error = np.linalg.norm(f - f_recovered_real)
    
    print(f"原始信號 f:\n {f}")
    print(f"逆變換 f' 的實部:\n {f_recovered_real.round(decimals=10)}")
    print(f"\n原始信號與逆變換結果的 L2 誤差: {error:.2e}")
    
    if error < 1e-9:
        print("結論: 誤差極小，驗證成功，DFT 和 IDFT 互為逆運算。")
    else:
        print("結論: 誤差較大，驗證失敗。")

3.2. 測試案例與執行
我們創建一個由兩個頻率分量組成的合成信號進行測試。

# 構造測試信號: 包含 5 Hz 和 12 Hz 兩個頻率分量
N_samples = 64
t = np.linspace(0, 1, N_samples, endpoint=False)
f_test = 2 * np.sin(2 * np.pi * 5 * t) + 1.5 * np.cos(2 * np.pi * 12 * t)

# 執行驗證
verify_dft_inverse(f_test)

3.3. 預期結果由於浮點數計算的極微小誤差，error 應該是一個非常小的數值 (通常在 $10^{-14}$ 範圍內)。

原始信號 $f$: 包含 64 個採樣點的合成正弦波\
恢復信號 $f'_{real}$: 數值上與 $f$ 完全一致 (四捨五入到小數點後 10 位)\
L2 誤差: 數值約為 $1.89\text{e-}14$ (或相似的極小值)\
結論: 誤差極小，驗證成功，DFT 和 IDFT 互為逆運算。
