# 離散傅立葉變換 (DFT) 與逆變換 (IDFT) 實現及可逆性驗證

這個專案使用 Python 和 NumPy 庫，從頭實現了離散傅立葉變換 (DFT) 及其逆變換 (IDFT)。核心目的是透過矩陣運算來高效計算這兩種變換，並驗證它們互為逆運算 (即 $IDFT(DFT(x)) \approx x$) 的數學性質。

## 1. 數學原理

### A. 離散傅立葉變換 (DFT)

DFT 將一個長度為 $N$ 的時域離散信號 $x[n]$ 轉換為頻域信號 $X[k]$。

$$X[k] = \sum_{n=0}^{N-1} x[n] \cdot e^{-i 2 \pi k n / N}$$

### B. 離散傅立葉逆變換 (IDFT)

IDFT 將頻域信號 $X[k]$ 逆轉換回原始時域信號 $x[n]$。

$$x[n] = \frac{1}{N} \sum_{k=0}^{N-1} X[k] \cdot e^{+i 2 \pi k n / N}$$

## 2. 程式碼實現 (使用 NumPy 矩陣運算)

本實現避免了多重迴圈，利用 NumPy 的向量化操作和矩陣乘法 (`.dot()`) 來計算變換矩陣 $M$ 與信號向量的乘積。

### 2.1. DFT 函數 (`dft`)

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

def idft(X):
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
