# 📐 Geometry Toolkit (Python) - 二維幾何計算工具庫

本專案實作了一組**高效且穩定的二維幾何計算工具庫**，涵蓋點、直線、圓等基本幾何物件的表示方式，以及求解常見幾何問題的演算法。

內容適合作為幾何計算、電腦圖學、競賽程式設計（競賽幾何）等用途的基礎工具。

## 📌 功能總覽 (Features)

### 幾何物件 (Geometry Objects)

| 類別 (Class) | 描述 (Description) | 表示法 (Representation) |
| :--- | :--- | :--- |
| `Point(x, y)` | 點 | $P(x, y)$ |
| `Line(a, b, c)` | 一般式直線 | $ax + by + c = 0$ |
| `Circle(center, r)` | 圓 | 中心點 `center`，半徑 `r` |
| `Triangle(A, B, C)` | 三角形 | 三個頂點 $A, B, C$ |

### 幾何運算 (Geometric Operations)

* `line_intersection(L1, L2)`: **兩直線交點**
* `line_circle_intersection(line, circle)`: **直線與圓交點**
* `circle_intersection(c1, c2)`: **兩圓交點**
* `foot_of_perpendicular(line, p)`: **一點到直線垂足** (投影點)
* `verify_pythagoras(A, B, C)`: **畢氏定理驗證** (檢查 $\angle C$ 是否為直角)

### 幾何變換 (Geometric Transformations)

所有變換函數都以新的 `Point` 物件回傳結果。

* `translate_point(p, dx, dy)`: **平移**
* `scale_point(p, sx, sy)`: **縮放**
* `rotate_point(p, theta)`: **旋轉**（以原點 $(0, 0)$ 為中心，角度 $\theta$ 以弧度 (radian) 計算）

---

## 🧩 解題流程與設計理念 (Design Philosophy)

本工具庫遵循由簡到繁的原則，以穩健的幾何和代數方法構築計算框架。

### 1. 定義基本幾何物件

* 以 Python 類別 (`Point`, `Line`, `Circle`) 封裝數據和行為，使操作更具物件導向性。
* 例如：`Line.from_points(p1, p2)` 可以從兩個點自動計算出直線的一般式 $ax + by + c = 0$。

### 2. 兩直線交點

* **方法：** 採用二元一次方程組的標準解法（如克拉瑪公式）。
    $$
    \begin{cases}
    a_1x + b_1y + c_1 = 0 \\
    a_2x + b_2y + c_2 = 0
    \end{cases}
    $$
* **判斷：** 計算判別式 $D = a_1b_2 - a_2b_1$。若 $D=0$，表示兩直線平行或重合，無單一交點。

### 3. 直線與圓交點

* **方法：** 將直線方程代入圓方程，化為**一元二次方程** $Ax^2 + Bx + C = 0$，再利用判別式 $\Delta = B^2 - 4AC$ 求解。
* **判別：**
    * $\Delta < 0$：無交點 (分離)
    * $\Delta = 0$：相切 (單一交點)
    * $\Delta > 0$：相交 (兩交點)
* **優化：** 實作上採用代數法結合幾何方式簡化（將圓心平移至原點），以提高計算穩定性。

### 4. 兩圓交點

採用**經典幾何推導法**，計算穩定且避免複雜的代數運算：

1.  計算兩圓心距離 $d$。
2.  處理特殊情況：外離、內含、同心，這些情況無交點。
3.  找到連心線上的基準點 $M$。
4.  利用 $M$ 點和垂直於連心線的方向向量，推算並得出兩交點座標。

### 5. 垂足（投影）計算

* **方法：** 針對點 $P(x_0, y_0)$ 投影到直線 $ax + by + c = 0$。
* **原理：** 公式源自向量投影以及直線的法向量 $\vec{n} = (a, b)$ 的性質。

### 6. 畢氏定理驗證

* **方法：** 驗證 $AB^2 = AC^2 + BC^2$。
* **關鍵：** 運算中**避免開根號** (即使用距離平方)，以最小化浮點數運算帶來的精度誤差。若等式成立，則 $\angle C$ 為直角。

### 7. 幾何變換（旋轉）

* **方法：** 採用標準的二維旋轉矩陣。
    $$
    \begin{pmatrix}
    x' \\
    y'
    \end{pmatrix}
    =
    \begin{pmatrix}
    \cos\theta & -\sin\theta \\
    \sin\theta & \cos\theta
    \end{pmatrix}
    \begin{pmatrix}
    x \\
    y
    \end{pmatrix}
    $$

## 🛠️ 使用環境

* **語言：** Python 3.x
* **依賴：** 僅依賴標準數學函式庫 `math` (用於三角函數和開根號)。
