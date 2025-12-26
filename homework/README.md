# 習題 1 : 請用程式驗證微積分基本定理

[習題](https://github.com/wengziting00/_cm/blob/main/homework/hw1/hw1.py)
[說明](https://github.com/wengziting00/_cm/blob/main/homework/hw1/README.md)
AI資料遺失

設定一個非常小的數 ℎ=0.00001 h=0.00001，用來模擬微分與積分的極小變化量。接著利用「差分法」定義 df(f, x)，以 (f(x+h)−f(x))/h 來近似函數在 x 處的導數；再用「黎曼和」的概念定義 integral(f, a, b)， 從 a 開始以步長 h 累加每一小段的面積 f(x)⋅h，近似計算定積分。之後在 theorem1(f, x) 中，先把函數 f 從 0 積分到 x，再對這個積分結果做微分，理論上應該回到原函數f(x)， 這正是在數值方式下驗證微積分基本定理。程式最後以 𝑓(𝑥)=x^3為例，計算在 x=2 處的導數近似值與從 0 到 2 的定積分近似值，並檢查「先積分再微分」的結果是否與原函數在該點的值足夠接近（誤差小於 0.01）

# 習題 2 : 請寫程式求解二次多項式的根
[習題二](https://github.com/wengziting00/_cm/blob/main/homework/hw2/hw2.py)
[說明](https://github.com/wengziting00/_cm/blob/main/homework/hw2/README.md)
AI資料遺失

首先匯入 cmath，因為當判別式 𝑏^2−4𝑎𝑐b2為負數時，平方根會是複數，cmath 可以正確處理這種情況。接著在 root2(a, b, c) 函式中，先計算判別式的平方根 discriminant = cmath.sqrt(b**2 - 4ac)，再依照二次方程式公式 2a分之−𝑏±根號b^2−4𝑎𝑐分別求出兩個解 root1 和 root2。為了確認計算是否正確，程式把這兩個解再代回原本的函數 𝑓(𝑥)=𝑎𝑥2+𝑏𝑥+𝑐 得到 f1 和 f2，理論上結果應該非常接近 0。最後使用 cmath.isclose 檢查代入後的結果是否在容許誤差範圍內，如果不接近 0 就印出警告訊息，提醒可能有數值誤差或計算問題，確認無誤後再把兩個根回傳

# 習題 3 : 請寫程式求解三次多項式的根 (加分題）
[習題三](https://github.com/wengziting00/_cm/blob/main/homework/hw3/hw3.py)
[說明](https://github.com/wengziting00/_cm/blob/main/homework/hw3/README.md)
AI資料遺失

先檢查 a!=0，再把方程式除以 (a) 使首項係數變成 1，接著透過代換 x=t−3分之b​ 將原式化為沒有二次項的簡化形式 t^3+pt+q=0，並由此計算 (p)、(q) 與判別式 Δ。之後套用 Cardano 公式，計算 u^3=-2分之q+根號Δ 與 u^3=-2分之q-根號Δ，再用複數立方根與立方根單位根 𝜔 ω 組合出三個 t 的解，最後將結果回代 x=t-3分之b，得到原三次方程式的三個解（包含實根與複根）

# 習題 4 : （思考）請寫一個函數 root(c) 求出 n 次多項式的根 （ n>=5 的時候，數學上證明沒有公式 -- 伽羅瓦定理）
[習題4](https://github.com/wengziting00/_cm/blob/main/homework/hw4/hw4.py)
[說明](https://github.com/wengziting00/_cm/blob/main/homework/hw4/README.md)
[AI](https://gemini.google.com/app/5295fd37c17d5673)

函式 solve_polynomial_roots(c) 接收一個串列 c，其中第 (i) 個元素代表 (x^i) 的係數，程式會先從最高次開始移除尾端多餘的 0，避免因為最高次係數為 0 而誤判多項式次數；如果移除後發現整個串列為空，表示多項式是 (f(x)=0)，此時任何 (x) 都是解，程式會印出警告並回傳空陣列。接著將係數順序反轉，因為 numpy.roots 需要由最高次到最低次的係數排列，然後印出原始係數與轉換後的係數供檢查。最後呼叫 np.roots 計算多項式的所有根（包含實根與複根），並將結果回傳。

# 第二週習題：有限體
[習題5](https://github.com/wengziting00/_cm/blob/main/homework/week2/week2.py)
[說明](https://github.com/wengziting00/_cm/blob/main/homework/week2/README.md)
[AI](https://gemini.google.com/app/513b61f539777e4d?hl=zh-TW)

## 有限體觀念簡介

有限體（Galois Field, GF）是一種代數結構，其運算規則與我們熟悉的數學運算相同，但元素個數有限。

### 體的三大公理

1. **加法形成交換群**（封閉、結合律、單位元、反元素、交換律）
2. **非零元素乘法形成交換群**
3. **乘法對加法滿足分配律**

### GF(5) 特性

- 元素：`{0, 1, 2, 3, 4}`
- 所有運算皆在 **模 5** 下進行（$ \bmod\ 5 $）
- 本專案透過 `GF5` 類別實作了：
  - 加法 `add(a, b)`
  - 乘法 `mul(a, b)`
  - 減法 `sub(a, b)`
  - 除法 `div(a, b)`

---

##  程式功能說明

### 檔案結構

| 檔案名稱            | 說明內容                                 |
|---------------------|------------------------------------------|
| `group_axioms.py`   | 驗證加法群與乘法群是否滿足群的五大公理     |
| `field_axioms.py`   | 驗證體的第三公理：乘法對加法的分配律       |
| `field_gf5.py`      | 定義 GF(5) 的基本運算與公理驗證流程       |
| `gf_element.py`     | 定義 `GFElement` 類，支援 `+ - * /` 運算  |

###  驗證內容對應

| 類別 / 函數            | 對應數學結構                 | 驗證內容                         |
|------------------------|------------------------------|----------------------------------|
| `GF5AddGroup`          | $(\mathbb{Z}_5, +)$          | 加法群公理                       |
| `GF5MulGroup`          | $(\mathbb{Z}_5 \setminus \{0\}, \cdot)$ | 乘法群公理（非零元素） |
| `check_group_axioms()` | 群公理                        | 封閉性、結合律、單位元、反元素、交換律 |
| `check_distributivity()` | 體公理                       | 驗證分配律 $a(b+c) = ab + ac$     |
| `GFElement` 類別       | 運算物件封裝                 | 支援直觀的 `+ - * /` 運算         |

---

#  執行結果範例

以下為部分運算與驗證結果（模 5）：

加法：3 + 4 = 2
乘法：3 * 4 = 2
反元素：3^-1 = 2，因為 3 * 2 = 1
除法：4 / 3 = 3
分配律驗證：3 * (4 + 2) = 3 = 3 * 4 + 3 * 2 → ✅

# 第三週習題：幾何學：（點，線，圓）世界的建構
[week3](https://github.com/wengziting00/_cm/blob/main/homework/week3/homework3.py)
[說明](https://github.com/wengziting00/_cm/tree/main/homework/week3)
[AI](https://gemini.google.com/app/7fb82a480baf3cb0?hl=zh-TW)

## Geometry Toolkit (Python) - 二維幾何計算工具庫

本專案實作了一組**高效且穩定的二維幾何計算工具庫**，涵蓋點、直線、圓等基本幾何物件的表示方式，以及求解常見幾何問題的演算法。

內容適合作為幾何計算、電腦圖學、競賽程式設計（競賽幾何）等用途的基礎工具。

## 功能總覽 (Features)

### 幾何物件 (Geometry Objects)

| 類別 (Class) | 描述 (Description) | 表示法 (Representation) |
| :--- | :--- | :--- |
| `Point(x, y)` | 點 | P(x, y) |
| `Line(a, b, c)` | 一般式直線 | ax + by + c = 0 |
| `Circle(center, r)` | 圓 | 中心點 `center`，半徑 `r` |
| `Triangle(A, B, C)` | 三角形 | 三個頂點 A, B, C |

### 幾何運算 (Geometric Operations)

* `line_intersection(L1, L2)`: **兩直線交點**
* `line_circle_intersection(line, circle)`: **直線與圓交點**
* `circle_intersection(c1, c2)`: **兩圓交點**
* `foot_of_perpendicular(line, p)`: **一點到直線垂足** (投影點)
* `verify_pythagoras(A, B, C)`: **畢氏定理驗證** (檢查 C 角是否為直角)

### 幾何變換 (Geometric Transformations)

所有變換函數都以新的 `Point` 物件回傳結果。

* `translate_point(p, dx, dy)`: **平移**
* `scale_point(p, sx, sy)`: **縮放**
* `rotate_point(p, theta)`: **旋轉**（以原點 (0, 0) 為中心，角度 $\theta$ 以弧度 (radian) 計算）

---

## 解題流程與設計理念 (Design Philosophy)

本工具庫遵循由簡到繁的原則，以穩健的幾何和代數方法構築計算框架。

### 1. 定義基本幾何物件

* 程式以 `Point`、`Line`、`Circle` 作為基底，使後續計算能以物件的方式操作，而非純數值推算。
* 例如：`Line.from_points(p1, p2)` 可自動求出直線一般式。

### 2. 兩直線交點

* **方法：** 利用二元一次方程組解法（克拉瑪公式）求交點。
    * 方程組：a1*x + b1*y + c1 = 0 和 a2*x + b2*y + c2 = 0
* **判斷：** 若判別式 D = a1*b2 - a2*b1 為 0，則表示兩直線平行或重合，無交點。

### 3. 直線與圓交點

* **方法：** 將直線代入圓方程，化為一元二次方程。
* **判別式：**
    * 判別式 < 0：無交點
    * 判別式 = 0：相切
    * 判別式 > 0：兩交點
* **優化：** 採用代數法結合幾何方式簡化（平移圓心）。

### 4. 兩圓交點

* **原理：** 根據兩圓交點的經典幾何推導：
    1. 計算兩圓心距離 d。
    2. 處理特殊情況：外離、內含、同心皆無交點。
    3. 找到連心線上的基準點 M。
    4. 利用垂線方向推算兩交點。
* **優勢：** 該方法計算穩定且避免不必要的代數複雜度。

### 5. 垂足（投影）計算

* **方法：** 對點 P 投影到直線 ax + by + c = 0。
* **原理：** 公式來自向量投影與直線法向量 (a, b) 的性質。

### 6. 畢氏定理驗證

* **方法：** 利用距離平方（避免開根號誤差）驗證：
    * **AB² = AC² + BC²**
* **結果：** 若成立則 C 為直角頂點。

### 7. 幾何變換（旋轉）

* **方法：** 採用標準的二維旋轉矩陣。
* **運算式：** * **新 x' = x * cos(theta) - y * sin(theta)**
    * **新 y' = x * sin(theta) + y * cos(theta)**
* **說明：** 這裡 `theta` 為旋轉角度（弧度）。

# 第八週習題：機率統計 - 檢定背後的數學原理
[AI](https://gemini.google.com/app/05721409f016f69e?hl=zh-TW)

# 第九周習題：資訊理論
[HW8](https://github.com/wengziting00/_cm/blob/main/homework/week9/homework9.py)
[說明](https://github.com/wengziting00/_cm/tree/main/homework/week9)
[AI](https://gemini.google.com/app/32e30bc46b73129c?hl=zh-TW)

## 1. 連續丟 10000 次正面的機率
從最基本的部分開始：公平銅板正面的機率是 0.5。
如果要連續 10000 次都正面，機率就是：
P = 0.5^10000
算出來是：
≈ 9.33 × 10^-301
這比想像中還誇張地小，用科學記號才看得見，不然直接印會變成 0。
主要結論：即使是 1/2，一連乘 10000 次後就幾乎趨近零。

## 2. 用 log 來計算 0.5^10000
直接算 p^n 容易 underflow，所以改用 log identity：
log(p^n) = n * log(p)
兩種方式：
先算 0.5^10000 再取 log
直接算 10000 * log(0.5)
算出來結果一樣。
之後再把 log 用 exp() 還原，也能得到跟第 1 題一樣的機率。
重點就是：用 log 能避免浮點數下限問題，計算更穩定。

## 3. 熵、交叉熵、KL 散度的計算
使用的分佈：
P = [0.5, 0.3, 0.2]   # 真實分佈
Q = [0.4, 0.4, 0.2]   # 模型預測

分別計算：
H(P)：平均要多少 bit 才能表達 P
H(P,Q)：如果資料照 P 產生，但你用 Q 去編碼，需要多少 bit
KL(P‖Q)：量化 P 和 Q 的差距
結果符合預期：
KL(P‖Q) = H(P, Q) - H(P)

也就是：
Q 越接近 P → 交叉熵越小
Q 越偏離 P → 交叉熵越大

## 4. 驗證交叉熵在 Q = P 時最小
換一組新的分佈：
P      = [0.1, 0.4, 0.3, 0.2]
Q_good = [0.1, 0.4, 0.3, 0.2]  # 完全一樣
Q_bad  = [0.4, 0.1, 0.3, 0.2]  # 排列錯了

計算：
H(P, P)      # 完美模型
H(P, Q_bad)  # 不好的模型

輸出結果滿足：
H(P, P) ≤ H(P, Q_bad)
這再次驗證「交叉熵最小值出現在 Q = P」。
這就是為什麼神經網路訓練會用 cross entropy。

## 5. 7-4 漢明碼：編碼與解碼
這裡使用：
G（產生矩陣）做編碼
H（校驗矩陣）做解碼與錯誤偵測
✦ 編碼
給定 4-bit data：
c = d @ G  (mod 2)
就會得到 7-bit codeword。

解碼
收到 r 之後：
syndrome = r @ H^T (mod 2)
syndrome = 000 → 沒錯
syndrome ≠ 000 → 代表錯誤位置（1～7）
之後把該位 flip 回來即可修正單錯誤。

✦ 測試結果：

無錯誤 → 正確取回原始 4 bits
在第 5 bit 製造錯誤 → syndrome 指出 5 → 修正後恢復原資料
符合漢明碼「能偵錯 1、能修正 1」的特性。
## 6. 夏農信道編碼定理 vs. 夏農–哈特利定理

這部分整理兩個定理的用途差異。
夏農信道編碼定理（Channel Coding Theorem）
重點是：
R < C → 有可能做到可靠傳輸（錯誤率 → 0）
R > C → 無論如何都不可能可靠傳輸
這是「到底能不能可靠通訊」的根本界線。

夏農–哈特利定理（Shannon–Hartley Theorem）
專門算 頻寬受限 + AWGN（高斯白噪聲）信道 的容量：
C = B * log2(1 + S/N)
容量由以下決定：
B（頻寬）SNR（訊噪比）
這是特定模型的「可達最高傳輸速度」。

# 第十周習題：線性代數
[HW9](https://github.com/wengziting00/_cm/blob/main/homework/week10/homework10.py)
[說明](https://github.com/wengziting00/_cm/blob/main/homework/week10/README.md)
[AI](https://gemini.google.com/app/5a46793b4b1d78a6?hl=zh-TW)

## 1. 什麼是「線性」？為何叫「代數」？

### 線性（Linearity）
一個運算 \(T\) 若滿足  
\[
T(au + bv) = aT(u) + bT(v)
\]  
則稱線性；它保持比例與可加性。

### 代數（Algebra）
線性代數研究向量、矩陣、線性變換，以及其運算規則，因此稱為「代數」。

---

## 2. 空間與向量空間

### 數學中的空間（Space）
空間＝元素集合＋定義其上的結構（如加法、距離、內積）。

### 向量空間（Vector Space）
能做「向量加法」與「數乘」，並滿足公理，因此被視為一種空間。

---

## 3. 矩陣與向量的關係

- 向量：座標表示。
- 矩陣：線性變換在特定基底下的座標表示。
- 矩陣的每欄＝基底向量經線性變換後的位置。

矩陣可表示旋轉、縮放、剪切、反射、投影等線性變換。

---

## 4. 矩陣如何表示 2D / 3D 的平移、縮放、旋轉？

### 2D 縮放

$$
\begin{bmatrix}
s_x & 0 \\
0 & s_y
\end{bmatrix}
$$

### 2D 旋轉矩陣

$$
\begin{bmatrix}
\cos\theta & -\sin\theta\\
\sin\theta & \cos\theta
\end{bmatrix}
$$

### 2D 平移矩陣（需齊次座標）

$$
\begin{bmatrix}
1 & 0 & t_x\\
0 & 1 & t_y\\
0 & 0 & 1
\end{bmatrix}
$$

---

## 5. 行列式的意義、遞迴公式、與體積關係

### 幾何意義

det(A) =  空間體積伸縮的倍數

- $|\det A| = 1$ ：體積不變  
- $\det A = 0$ ：空間被壓扁（不可逆）

### Laplace 展開（遞迴公式）

$\det(A) = \sum_{j=1}^n (-1)^{1+j} a_{1j} M_{1j}$

### 與體積的關係

Volume after $= |\det A| \times$ Volume before


---

## 6. 對角化如何快速計算行列式

若

$$
A = P D P^{-1}
$$

則

$$
\det(A) = \det(D) = \prod_i \lambda_i
$$

即行列式為所有特徵值的乘積。

---

## 7. LU 分解快速計算行列式

若

$$
A = LU
$$

則

$$
\det(A) = \prod_i u_{ii}
$$

因為 $L$ 的對角線多為 1。

---

## 8. 特徵值與特徵向量的意義與用途

### 意義

$$
A v = \lambda v
$$

特徵向量方向不變，特徵值代表伸縮倍數。

### 用途
- 對角化  
- 計算 $A^n$  
- Markov 鏈  
- PCA  
- 模式識別  
- 微分方程與物理振動


---

## 9. QR 分解

$$
A = QR
$$

- \(Q\)：正交矩陣  
- \(R\)：上三角矩陣  

用途：解方程、數值穩定、求特徵值的基礎。

---

## 10. 用 QR 分解求特徵值（QR Algorithm）

重複：

1. $A_k = Q_k R_k$
2. $A_{k+1} = R_k Q_k$

最終  $A_k$ 會收斂到上三角矩陣，其對角線即特徵值。

---

## 11. SVD 與特徵值分解的關係

### SVD
$$
A = U \Sigma V^T
$$

- \(U, V\)：正交  
- $\Sigma$：奇異值  

### 與特徵值關係

$$
A^T A = V \Sigma^2 V^T
$$



奇異值＝ $(A^T A)$ 的特徵值的平方根。

---

## 12. PCA 與 SVD 的關係

PCA：找出資料中變異最大的方向（降維）。

資料矩陣 \(X\)：

$X = U \Sigma V^T$

- \(V\)：主成分方向  
- \(\Sigma\)：變異大小（奇異值）  
- 選取前 k 個奇異值 → k 維 PCA  

---

# 第11周習題：請寫出傅立葉正轉換和逆轉換的函數（不要用套件）

[習題10](https://github.com/wengziting00/_cm/blob/main/homework/week11/HW10.py)
[說明](https://github.com/wengziting00/_cm/tree/main/homework/week11)
[AI](https://gemini.google.com/app/5a46793b4b1d78a6?hl=zh-TW)
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
```
 2.2. IDFT 函數 (idft)在 IDFT 中，變換矩陣的指數符號相反 (變為正號)，並且結果需要除以信號長度 $N$。

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

## 3. 可逆性驗證verify_dft_inverse 函數執行 $DFT \to IDFT$ 的過程，並計算原始信號與恢復信號之間的 L2 誤差。

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
N_samples = 64\
t = np.linspace(0, 1, N_samples, endpoint=False)\
f_test = 2 * np.sin(2 * np.pi * 5 * t) + 1.5 * np.cos(2 * np.pi * 12 * t)

# 執行驗證
verify_dft_inverse(f_test)

3.3. 預期結果由於浮點數計算的極微小誤差，error 應該是一個非常小的數值 (通常在 $10^{-14}$ 範圍內)。

| 項目 | 描述 |
| :--- | :--- |
| 原始信號 f | 包含 64 個採樣點的合成正弦波 |
| 恢復信號 f'real | 數值上與 f 完全一致 (四捨五入到小數點後 10 位) |
| L2 誤差 | 數值約為 1.89e-14 (或相似的極小值) |
| 結論 | 誤差極小，驗證成功 |

# 第13周習題：請寫程式求解常係數齊次常微分方程
[習題11](https://github.com/wengziting00/_cm/blob/main/homework/week13/week13.py)
[說明](https://github.com/wengziting00/_cm/blob/main/homework/week13/README.md)
[AI](https://gemini.google.com/app/c5f67c12dd438c99?hl=zh-TW)

根據微分方程理論，這類方程的通解完全由其特徵方程的根所決定，因此程式首先利用 numpy.roots 計算由輸入係數所形成的特徵多項式的所有根。由於數值計算可能產生微小誤差，例如理論上的實根在計算後可能帶有極小的虛部，程式會對根進行修正，將虛部接近零的根視為實根，以確保分類的正確性。

接著，程式利用 Counter 統計各特徵根的重複次數，因為在常係數齊次線性微分方程中，重根的次數會直接影響通解中所需乘上的多項式因子。對於實根，程式依據重根次數產生 𝑒𝜆𝑥,𝑥𝑒𝜆𝑥,𝑥2𝑒𝜆𝑥eλx,xeλx,x2eλx 等解；對於複數共軛根，則依理論轉換為以 𝑒𝛼𝑥cos⁡(𝛽𝑥)αxcos(βx) 與 𝑒𝛼𝑥sin⁡(𝛽𝑥)eαxsin(βx) 為基底的實值解。如果複數根為重根，則同樣乘上對應次數的 𝑥𝑘xk

最後，程式將所有線性獨立解以常數 𝐶𝑖Ci組合成通解，並以字串形式輸出
