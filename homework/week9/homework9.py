# 1. 計算公平銅板連續投擲 10000 次，全部得到正面的機率 (p^10000)

p = 0.5
n = 10000

probability = p**n
print(f"公平銅板連續投擲 {n} 次，全部得到正面的機率 (0.5^{n}):")
# 使用科學記號顯示以避免顯示為 0.0
print(f"P = {probability:.4e}")
# 輸出：P = 9.3326e-301

import numpy as np


# 2. 使用 log(p^n) = n log(p) 計算 log(0.5^10000)

p = 0.5
n = 10000

# 方法一：先計算 p^n, 再取 log
log_pn_method1 = np.log(p**n)

# 方法二：使用 n * log(p)
log_p = np.log(p) # np.log() 預設為自然對數 ln
log_pn_method2 = n * log_p

# 為了驗證結果，我們可以將 log(p^n) 轉換回機率 P = e^(log(p^n))
probability_check = np.exp(log_pn_method2)

print(f"\n計算 log(0.5^{n}):")
print(f"n * log(p) = {n} * {log_p:.4f} = {log_pn_method2:.4f}")
print(f"log(p^n) (方法一) = {log_pn_method1:.4f}")

print(f"\n從 log 值還原機率: exp({log_pn_method2:.4f}) = {probability_check:.4e} (與第1題結果一致)")

import numpy as np

# 3. 計算『熵，交叉熵，KL 散度』
# 為了避免 log(0)，我們在機率為 0 時忽略該項，或加上一個極小值 (EPSILON)

def entropy(P):
    """計算機率分佈 P 的熵 H(P) (以 2 為底)"""
    P = np.asarray(P)
    # 忽略 P=0 的情況，因為 p * log(p) -> 0
    P = P[P > 1e-12]
    return -np.sum(P * np.log2(P))

def cross_entropy(P, Q):
    """計算機率分佈 P 和 Q 的交叉熵 H(P, Q) (以 2 為底)"""
    P = np.asarray(P)
    Q = np.asarray(Q)
    # 忽略 P=0 的情況
    mask = P > 1e-12
    return -np.sum(P[mask] * np.log2(Q[mask]))

def kl_divergence(P, Q):
    """計算 KL 散度 D_KL(P || Q) (以 2 為底)"""
    P = np.asarray(P)
    Q = np.asarray(Q)
    # 忽略 P=0 的情況
    mask = P > 1e-12
    # Q[mask] 不能為 0
    Q_safe = Q[mask]
    if np.any(Q_safe < 1e-12):
        print("警告: D_KL(P || Q) 無窮大 (Q 中的某些項為 0)")
        return np.inf

    return np.sum(P[mask] * np.log2(P[mask] / Q_safe))


# 範例機率分佈 (假設是一個三種結果的事件)
P_true = np.array([0.5, 0.3, 0.2]) # 真實分佈 (True distribution)
Q_model = np.array([0.4, 0.4, 0.2]) # 模型預測分佈 (Model distribution)

print("\n--- 3. 資訊理論量計算 ---")
H_P = entropy(P_true)
H_PQ = cross_entropy(P_true, Q_model)
D_KL_PQ = kl_divergence(P_true, Q_model)

print(f"真實分佈 P: {P_true}")
print(f"預測分佈 Q: {Q_model}")
print(f"1. 熵 H(P) (資訊量的期望): {H_P:.4f} bits")
print(f"2. 交叉熵 H(P, Q) (模型預測 P 的平均位元數): {H_PQ:.4f} bits")
print(f"3. KL 散度 D_KL(P || Q) (兩個分佈的差異): {D_KL_PQ:.4f} bits")
# 驗證 KL = Cross_Entropy - Entropy
print(f"   驗證: H(P, Q) - H(P) = {H_PQ - H_P:.4f} (與 D_KL 一致)")

# 4. 驗證 cross_entropy(P, P) <= cross_entropy(P, Q)

P = np.array([0.1, 0.4, 0.3, 0.2]) # 真實分佈 P (已排序，但總和為 1)
Q_good = np.array([0.1, 0.4, 0.3, 0.2]) # 完美的模型 Q = P
Q_bad = np.array([0.4, 0.1, 0.3, 0.2]) # 較差的模型 Q != P

# 計算 H(P, P)
H_PP = cross_entropy(P, Q_good)
# 計算 H(P, Q)
H_PQ = cross_entropy(P, Q_bad)

print("\n--- 4. 交叉熵比較驗證 ---")
print(f"分佈 P: {P}")
print(f"分佈 Q (完美): {Q_good}")
print(f"分佈 Q (較差): {Q_bad}")

print(f"H(P, P) (最小交叉熵/熵): {H_PP:.4f}")
print(f"H(P, Q) (較差交叉熵):    {H_PQ:.4f}")

if H_PP <= H_PQ:
    print(f"驗證結果: {H_PP:.4f} <= {H_PQ:.4f}，**驗證成功** (交叉熵在 Q=P 時最小)")
else:
    print(f"驗證結果: {H_PP:.4f} > {H_PQ:.4f}，**驗證失敗**")

import numpy as np

# 5. 『7-4 漢明碼』的編碼與解碼程式 (使用 Modulo 2 運算)

# 產生矩陣 (Generator Matrix G)
G = np.array([
    [1, 0, 0, 0, 1, 1, 0],
    [0, 1, 0, 0, 1, 0, 1],
    [0, 0, 1, 0, 0, 1, 1],
    [0, 0, 0, 1, 1, 1, 1]
])

# 校驗矩陣 (Parity Check Matrix H) - 用於解碼
H = np.array([
    [1, 1, 0, 1, 1, 0, 0],
    [1, 0, 1, 1, 0, 1, 0],
    [0, 1, 1, 1, 0, 0, 1]
])


def encode_hamming_7_4(data_word):
    """7-4 漢明碼編碼: c = d * G mod 2"""
    data_word = np.array(data_word)
    code_word = (data_word @ G) % 2
    return code_word

def decode_hamming_7_4(received_word):
    """7-4 漢明碼解碼：計算伴隨式 s = r * H^T mod 2"""
    received_word = np.array(received_word)
    # 計算伴隨式 (Syndrome)
    # H.T 是 H 的轉置
    syndrome = (received_word @ H.T) % 2

    syndrome_int = int("".join(map(str, syndrome)), 2)

    # 伴隨式對應到錯誤位置 (H 的列向量剛好是 1 到 7 的二進位表示)
    if syndrome_int != 0:
        error_pos = syndrome_int
        print(f"  偵測到錯誤！伴隨式: {syndrome} (十進位 {error_pos})")
        
        # 校正錯誤 (將該位置取反 0->1, 1->0)
        corrected_word = received_word.copy()
        # 錯誤位置從 1 開始，陣列索引從 0 開始
        corrected_word[error_pos - 1] = 1 - corrected_word[error_pos - 1]
        
        # 提取資料位元 (碼字的前 4 位元)
        decoded_data = corrected_word[:4]
        return decoded_word, True # 校正後的碼字, 已校正標誌
    else:
        # 無錯誤，直接提取資料位元
        decoded_data = received_word[:4]
        return decoded_data, False # 解碼資料, 未校正標誌

# --- 範例執行 ---
data = [1, 0, 1, 1] # 4 位元資料
code_word = encode_hamming_7_4(data)
print("\n--- 5. 7-4 漢明碼 ---")
print(f"原始資料 (d): {data}")
print(f"編碼結果 (c): {code_word}")

# 1. 無錯誤傳輸
received_no_error = code_word.copy()
decoded_no_error, _ = decode_hamming_7_4(received_no_error)
print(f"\n[無錯誤] 接收碼字: {received_no_error} -> 解碼資料: {decoded_no_error}")


# 2. 引入一個錯誤 (第 5 個位元錯誤)
error_pos = 5
received_error = code_word.copy()
received_error[error_pos - 1] = 1 - received_error[error_pos - 1] # 反轉第 5 位
print(f"\n[單錯誤] 接收碼字: {received_error}")
decoded_error, is_corrected = decode_hamming_7_4(received_error)
print(f"  校正後資料: {decoded_error} (與原始資料 {data} 比對: {np.array_equal(decoded_error, data)})")

"""6. 說明『夏農信道編碼定理』和『夏農-哈特利定理 (Shannon–Hartley Theorem)』
夏農信道編碼定理（Channel Coding Theorem）
告訴我們：只要傳輸速率 R 小於信道容量 C，就能用適當編碼讓錯誤率趨近 0。
若 R > C，再怎麼編碼都無法可靠傳輸。
是「可靠通訊是否可能」的普遍原則。

夏農–哈特利定理（Shannon–Hartley Theorem）
專門針對 頻寬受限、AWGN 高斯噪聲信道 的容量計算公式：C=Blog2​(1+S/N)
告訴我們容量受 頻寬 B 與 SNR 決定。"""
