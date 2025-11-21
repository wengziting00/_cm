# 📐 Geometry Toolkit (Python) - 二維幾何計算工具庫

本專案實作了一組**高效且穩定的二維幾何計算工具庫**，涵蓋點、直線、圓等基本幾何物件的表示方式，以及求解常見幾何問題的演算法。

內容適合作為幾何計算、電腦圖學、競賽程式設計（競賽幾何）等用途的基礎工具。

## 📌 功能總覽 (Features)

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

## 🧩 解題流程與設計理念 (Design Philosophy)

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
* **優化：** 採用代數法結合幾何方式簡
