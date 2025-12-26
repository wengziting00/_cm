# finite_field_gf5_simple.py

# =========================================================
# 參數設定：使用 GF(5)，即質數 P = 5
# =========================================================
P = 5

class GFElement:
    """
    有限體 GF(P) 的元素類別 (P=5)，支援運算子重載。
    這個類別同時是我們進行群公理和體公理驗證的基礎。
    """
    
    def __init__(self, value):
        if not isinstance(value, int):
             raise TypeError(f"GFElement 必須是整數")
        self._value = value % P

    @property
    def value(self):
        return self._value
        
    def __repr__(self):
        return f"GF{P}({self.value})"

    def __eq__(self, other):
        """判斷相等，支援與整數比較"""
        if isinstance(other, GFElement):
            return self.value == other.value
        return self.value == (other % P)
        
    def __hash__(self):
        return hash(self.value)
        
    # --- 運算子重載 ---

    def _inverse_mul(self):
        """計算乘法反元素 (僅供內部使用)"""
        if self.value == 0:
            raise ZeroDivisionError("零沒有乘法反元素")
        # 尋找 x 使得 self.value * x = 1 (mod P)
        for x in range(1, P):
            if (self.value * x) % P == 1:
                return x
        # 對於質數 P，這行不應該被執行
        raise RuntimeError(f"無法找到 {self.value} 的乘法反元素")

    def __add__(self, other):
        other_val = other.value if isinstance(other, GFElement) else other
        return GFElement((self.value + other_val) % P)

    def __sub__(self, other):
        other_val = other.value if isinstance(other, GFElement) else other
        return GFElement((self.value - other_val) % P)

    def __mul__(self, other):
        other_val = other.value if isinstance(other, GFElement) else other
        return GFElement((self.value * other_val) % P)

    def __truediv__(self, other):
        other_val = other.value if isinstance(other, GFElement) else other
        
        if other_val == 0:
            raise ZeroDivisionError("除以零")
            
        # 進行除法： self / other = self * other^-1
        other_inv = GFElement(other_val)._inverse_mul()
        return GFElement(self.value * other_inv)

    # 支援右側運算 (例如 2 + GF5(3))
    __radd__ = __add__
    __rmul__ = __mul__


# =========================================================
# 抽象群結構 (用於公理驗證)
# 這裡我們使用 GFElement 類別的運算來定義群操作
# =========================================================

# 所有的元素集合 {0, 1, 2, 3, 4}
ALL_ELEMENTS = [GFElement(i) for i in range(P)]
# 非零元素集合 {1, 2, 3, 4}
NON_ZERO_ELEMENTS = [GFElement(i) for i in range(1, P)]

class GFAddGroup:
    """GF(P) 的加法群 (F, +) 結構"""
    @staticmethod
    def op(a, b): return a + b
    @staticmethod
    def identity(): return GFElement(0)
    @staticmethod
    def inverse(a): return GFElement(0) - a
        
class GFMulGroup:
    """GF(P) 非零元素的乘法群 (F \ {0}, *) 結構"""
    @staticmethod
    def op(a, b): return a * b
    @staticmethod
    def identity(): return GFElement(1)
    @staticmethod
    def inverse(a): return GFElement(1) / a # 利用運算子重載的 / 來計算反元素


# =========================================================
# 公理驗證函數
# =========================================================

def check_group_axioms(GroupClass, elements, name):
    """驗證群公理 (封閉性、結合律、單位元素、反元素、交換律)"""
    
    identity = GroupClass.identity()
    element_set = set(elements)
    
    print(f"\n--- 驗證 {name} (運算: {GroupClass.op(GFElement(1), GFElement(2))}) ---")
    
    # 1. 封閉性
    for a in elements:
        for b in elements:
            c = GroupClass.op(a, b)
            if c not in element_set:
                print(f"失敗: 封閉性不成立: {a} op {b} = {c}, 但 {c} 不在集合中。")
                return False
    print("✓ 封閉性成立")

    # 2. 結合律
    for a in elements:
        for b in elements:
            for c in elements:
                left = GroupClass.op(GroupClass.op(a, b), c)
                right = GroupClass.op(a, GroupClass.op(b, c))
                if left != right:
                    print(f"失敗: 結合律不成立: ({a} op {b}) op {c} = {left}, 而 {a} op ({b} op {c}) = {right}")
                    return False
    print("✓ 結合律成立")

    # 3. 單位元素
    for a in elements:
        if GroupClass.op(a, identity) != a or GroupClass.op(identity, a) != a:
            print(f"失敗: 單位元素 {identity} 不成立: {a} op {identity} != {a}")
            return False
    print(f"✓ 單位元素 {identity} 存在")

    # 4. 反元素
    for a in elements:
        inverse = GroupClass.inverse(a)
        if inverse not in element_set:
            print(f"失敗: 元素 {a} 的反元素 {inverse} 不在集合中。")
            return False
        if GroupClass.op(a, inverse) != identity or GroupClass.op(inverse, a) != identity:
            print(f"失敗: 反元素不成立: {a} op {inverse} != {identity}")
            return False
    print("✓ 反元素存在")
    
    # 5. 交換律
    for a in elements:
        for b in elements:
            if GroupClass.op(a, b) != GroupClass.op(b, a):
                print(f"失敗: 交換律不成立: {a} op {b} != {b} op {a}")
                return False
    print("✓ 交換律成立 (交換群)")

    print(f"*** {name} 是一個交換群 (通過) ***")
    return True

def check_distributivity(elements):
    """驗證乘法對加法的分配律"""
    
    print("\n--- 驗證 分配律 (Field Axiom) ---")
    
    for a in elements:
        for b in elements:
            for c in elements:
                # a * (b + c) = a * b + a * c
                left_side = a * (b + c)
                right_side = (a * b) + (a * c)
                
                if left_side != right_side:
                    print(f"失敗: 分配律不成立: {a} * ({b} + {c}) = {left_side}, 而 {a} * {b} + {a} * {c} = {right_side}")
                    return False
                    
    print("✓ 分配律成立: a * (b + c) = a * b + a * c")
    print("*** 體公理驗證通過 (分配律) ***")
    return True


# =========================================================
# 主程式：執行驗證與示範
# =========================================================

def main():
    print(f"==========================================")
    print(f"  有限體 GF({P}) 公理驗證與示範 (單一檔案)  ")
    print(f"==========================================")
    
    # 1. 群公理驗證
    check_group_axioms(GFAddGroup, ALL_ELEMENTS, f"GF({P}) 加法群")
    check_group_axioms(GFMulGroup, NON_ZERO_ELEMENTS, f"GF({P}) 非零乘法群")
    
    # 2. 體公理 (分配律) 驗證
    check_distributivity(ALL_ELEMENTS)
    
    # 3. 運算子重載示範 (類似整數操作)
    print(f"\n==========================================")
    print(f"  GF({P}) 運算子重載示範 ")
    print(f"==========================================")

    a = GFElement(3)
    b = GFElement(4)
    c = GFElement(2)

    # 加法: 3 + 4 = 7 = 2 (mod 5)
    d_add = a + b 
    print(f"加法: {a} + {b} = {d_add}")

    # 乘法: 3 * 4 = 12 = 2 (mod 5)
    d_mul = a * b
    print(f"乘法: {a} * {b} = {d_mul}")
    
    # 減法: 2 - 4 = -2 = 3 (mod 5)
    d_sub = c - b
    print(f"減法: {c} - {b} = {d_sub}")
    
    # 除法: 4 / 3 = 4 * 3^-1 = 4 * 2 = 8 = 3 (mod 5)
    d_div = b / a
    print(f"除法: {b} / {a} = {d_div}")
    
    # 分配律示範: 3 * (4 + 2) = 3 * 4 + 3 * 2
    left = a * (b + c)
    right = (a * b) + (a * c)
    print(f"分配律示範: {a} * ({b} + {c}) = {left}")
    print(f"            ({a} * {b}) + ({a} * {c}) = {right}")
    print(f"結果確認: {left == right}")

if __name__ == "__main__":
    main()
