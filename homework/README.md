## 習題 1 : 請用程式驗證微積分基本定理

[習題](https://github.com/wengziting00/_cm/blob/main/homework/hw1/hw1.py)
[說明](https://github.com/wengziting00/_cm/blob/main/homework/hw1/README.md)
AI資料遺失

設定一個非常小的數 ℎ=0.00001 h=0.00001，用來模擬微分與積分的極小變化量。接著利用「差分法」定義 df(f, x)，以 (f(x+h)−f(x))/h 來近似函數在 x 處的導數；再用「黎曼和」的概念定義 integral(f, a, b)， 從 a 開始以步長 h 累加每一小段的面積 f(x)⋅h，近似計算定積分。之後在 theorem1(f, x) 中，先把函數 f 從 0 積分到 x，再對這個積分結果做微分，理論上應該回到原函數f(x)， 這正是在數值方式下驗證微積分基本定理。程式最後以 𝑓(𝑥)=x^3為例，計算在 x=2 處的導數近似值與從 0 到 2 的定積分近似值，並檢查「先積分再微分」的結果是否與原函數在該點的值足夠接近（誤差小於 0.01）

## 習題 2 : 請寫程式求解二次多項式的根
[習題二](https://github.com/wengziting00/_cm/blob/main/homework/hw2/hw2.py)
[說明](https://github.com/wengziting00/_cm/blob/main/homework/hw2/README.md)
AI資料遺失

首先匯入 cmath，因為當判別式 𝑏^2−4𝑎𝑐b2為負數時，平方根會是複數，cmath 可以正確處理這種情況。接著在 root2(a, b, c) 函式中，先計算判別式的平方根 discriminant = cmath.sqrt(b**2 - 4ac)，再依照二次方程式公式 2a分之−𝑏±根號b^2−4𝑎𝑐分別求出兩個解 root1 和 root2。為了確認計算是否正確，程式把這兩個解再代回原本的函數 𝑓(𝑥)=𝑎𝑥2+𝑏𝑥+𝑐 得到 f1 和 f2，理論上結果應該非常接近 0。最後使用 cmath.isclose 檢查代入後的結果是否在容許誤差範圍內，如果不接近 0 就印出警告訊息，提醒可能有數值誤差或計算問題，確認無誤後再把兩個根回傳

## 習題 3 : 請寫程式求解三次多項式的根 (加分題）
[習題三](https://github.com/wengziting00/_cm/blob/main/homework/hw3/hw3.py)
[說明](https://github.com/wengziting00/_cm/blob/main/homework/hw3/README.md)
AI資料遺失

先檢查 a!=0，再把方程式除以 (a) 使首項係數變成 1，接著透過代換 x=t−3分之b​ 將原式化為沒有二次項的簡化形式 t^3+pt+q=0，並由此計算 (p)、(q) 與判別式 Δ。之後套用 Cardano 公式，計算 u^3=-2分之q+根號Δ 與 u^3=-2分之q-根號Δ，再用複數立方根與立方根單位根 𝜔 ω 組合出三個 t 的解，最後將結果回代 x=t-3分之b，得到原三次方程式的三個解（包含實根與複根）

## 習題 4 : （思考）請寫一個函數 root(c) 求出 n 次多項式的根 （ n>=5 的時候，數學上證明沒有公式 -- 伽羅瓦定理）
[習題4](https://github.com/wengziting00/_cm/blob/main/homework/hw4/hw4.py)
[說明](https://github.com/wengziting00/_cm/blob/main/homework/hw4/README.md)
[AI](https://gemini.google.com/app/5295fd37c17d5673)

函式 solve_polynomial_roots(c) 接收一個串列 c，其中第 (i) 個元素代表 (x^i) 的係數，程式會先從最高次開始移除尾端多餘的 0，避免因為最高次係數為 0 而誤判多項式次數；如果移除後發現整個串列為空，表示多項式是 (f(x)=0)，此時任何 (x) 都是解，程式會印出警告並回傳空陣列。接著將係數順序反轉，因為 numpy.roots 需要由最高次到最低次的係數排列，然後印出原始係數與轉換後的係數供檢查。最後呼叫 np.roots 計算多項式的所有根（包含實根與複根），並將結果回傳。

## 第二週習題：有限體
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

##  執行結果範例

以下為部分運算與驗證結果（模 5）：

加法：3 + 4 = 2
乘法：3 * 4 = 2
反元素：3^-1 = 2，因為 3 * 2 = 1
除法：4 / 3 = 3
分配律驗證：3 * (4 + 2) = 3 = 3 * 4 + 3 * 2 → ✅



