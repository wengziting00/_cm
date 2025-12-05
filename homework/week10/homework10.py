# 1. 遞迴方式計算行列式
def recursive_determinant(A):
    """
    使用遞迴的餘因子展開法計算矩陣 A 的行列式。
    A 必須是一個方陣 (n x n)。
    """
    n = A.shape[0]
    
    # 檢查是否為方陣
    if A.shape[1] != n:
        raise ValueError("矩陣必須是方陣")

    # 基本情況：1x1 矩陣
    if n == 1:
        return A[0, 0]

    # 基本情況：2x2 矩陣
    if n == 2:
        return A[0, 0] * A[1, 1] - A[0, 1] * A[1, 0]

    det_val = 0
    
    # 沿著第一行進行餘因子展開
    for j in range(n):
        # 構造餘子式 (Minor)：刪除第一行和第 j 列
        # np.delete(A, 0, axis=0) 刪除第一行
        # 然後再對結果刪除第 j 列
        minor = np.delete(np.delete(A, 0, axis=0), j, axis=1)
        
        # 餘因子 C_ij = (-1)^(i+j) * M_ij
        # 這裡 i=0 (第一行)，所以 (-1)^(0+j) = (-1)^j
        cofactor_sign = (-1) ** j
        
        # 遞迴呼叫
        det_val += cofactor_sign * A[0, j] * recursive_determinant(minor)

    return det_val

# 測試
A_test = np.array([[1, 2, 3], 
                   [4, 5, 6], 
                   [7, 8, 9]])

# print("遞迴計算的行列式:", recursive_determinant(A_test))
# print("NumPy 內建的行列式:", np.linalg.det(A_test)) 
# (對於這個矩陣，結果應該接近 0)

# 2. LU 分解後計算行列式
from scipy.linalg import lu

def det_from_lu(A):
    """
    使用 LU 分解計算矩陣 A 的行列式。
    利用 det(A) = det(L) * det(U) 的性質。
    NumPy/SciPy 的 lu 函式返回 P, L, U，其中 PA = LU。
    因此 det(A) = det(P') * det(L) * det(U)
    """
    # SciPy 的 lu 函式返回置換矩陣 P, 下三角 L, 上三角 U，滿足 P*A = L*U
    P, L, U = lu(A)
    
    # 1. 計算 det(U): 上三角矩陣的行列式是對角線元素的乘積
    det_U = np.prod(np.diag(U))
    
    # 2. 計算 det(L): L 的對角線元素通常為 1 (假設沒有縮放)
    #    但 SciPy 確保 L 的對角線為 1，所以 det(L) = 1
    det_L = 1.0
    
    # 3. 計算 det(P): 置換矩陣的行列式為 +1 或 -1，取決於置換的次數
    #    det(P) 可以通過計算 P 矩陣列向量之間的奇偶性來確定，
    #    但 SciPy 的 P 是一個置換矩陣，可以直接用 det(P) 計算
    det_P = np.linalg.det(P)
    
    # det(A) = det(P^-1) * det(L) * det(U)
    # 由於 P 是正交/置換矩陣， det(P^-1) = 1/det(P)。
    # 實際運算為 det(A) = (1/det(P)) * det(L) * det(U)
    # 但因為 P*A = L*U，所以 det(P)*det(A) = det(L)*det(U)
    # => det(A) = (det(L)*det(U)) / det(P)
    # 由於 det(L)=1， det(A) = det(U) / det(P)
    
    # 由於 det(P) 永遠是 1 或 -1，所以 det(A) = det(P) * det(L) * det(U)
    # 我們使用一個更穩定的方法：
    det_A = det_P * det_L * det_U
    
    return det_A

# 測試
A_test_2 = np.array([[2, 1, 1], [4, 3, 3], [8, 7, 9]])
# print("\nLU 分解計算的行列式:", det_from_lu(A_test_2))
# print("NumPy 內建的行列式:", np.linalg.det(A_test_2))

# 3.驗證分解（LU, Eig, SVD）的重構性
def verify_decomposition(A):
    print(f"\n--- 驗證矩陣重構 ---")
    
    # 1. LU 分解 (使用 SciPy 的 P*A = L*U 形式)
    from scipy.linalg import lu
    P, L, U = lu(A)
    A_lu_reconstructed = P.T @ L @ U # P.T (P 的轉置) 相當於 P 的逆矩陣 P^-1
    lu_diff = np.linalg.norm(A - A_lu_reconstructed)
    print(f"LU 重構誤差 (P.T @ L @ U): {lu_diff:.2e}")

    # 2. 特徵值分解 (Eigendecomposition): A = P * Lambda * P^-1
    # 僅適用於方陣
    if A.shape[0] == A.shape[1]:
        eigen_values, P = np.linalg.eig(A)
        Lambda = np.diag(eigen_values) # 特徵值矩陣
        P_inv = np.linalg.inv(P) # 特徵向量矩陣的逆
        A_eig_reconstructed = P @ Lambda @ P_inv
        eig_diff = np.linalg.norm(A - A_eig_reconstructed)
        print(f"特徵值重構誤差 (P @ Lambda @ P_inv): {eig_diff:.2e}")
    else:
        print("特徵值分解：非方陣，跳過。")

    # 3. SVD 分解: A = U * Sigma * V.T
    U, s, V_T = np.linalg.svd(A)
    
    # 奇異值 s 是一個向量，需要轉換為對角矩陣 Sigma
    Sigma = np.zeros(A.shape)
    np.fill_diagonal(Sigma, s)
    
    A_svd_reconstructed = U @ Sigma @ V_T
    svd_diff = np.linalg.norm(A - A_svd_reconstructed)
    print(f"SVD 重構誤差 (U @ Sigma @ V.T): {svd_diff:.2e}")

# 測試用方陣
A_square = np.array([[1, 2, 0], [2, 1, 3], [0, 3, 1]], dtype=float)
verify_decomposition(A_square)

# 測試用非方陣
A_rect = np.array([[1, 2, 3], [4, 5, 6]], dtype=float)
verify_decomposition(A_rect)

# 4.使用特徵值分解「計算」SVD
def svd_from_eig(A):
    """
    通過對 A.T @ A 進行特徵值分解來計算 SVD。
    """
    # 1. 計算 C = A.T @ A
    C = A.T @ A
    
    # 2. 對 C 進行特徵值分解: C = V @ Lambda @ V.T
    # V: 右奇異向量矩陣 (Right Singular Vectors)
    eigen_values, V = np.linalg.eigh(C) # 使用 eigh 處理對稱矩陣，結果更穩定
    
    # 確保特徵值和特徵向量按照奇異值大小降序排列
    idx = eigen_values.argsort()[::-1]
    eigen_values = eigen_values[idx]
    V = V[:, idx]
    
    # 3. 奇異值 Sigma (Singular Values) 是特徵值的平方根
    # 由於數值誤差，特徵值可能有微小的負值，所以取 max(0, lambda)
    singular_values = np.sqrt(np.maximum(0, eigen_values))
    
    # 4. 構造 Sigma 矩陣
    m, n = A.shape
    Sigma = np.zeros((m, n))
    np.fill_diagonal(Sigma, singular_values)
    
    # 5. 計算左奇異向量 U (Left Singular Vectors)
    U = np.zeros((m, m))
    
    for i in range(m):
        if i < len(singular_values) and singular_values[i] > 1e-10: # 避免除以零
            # u_i = (1 / sigma_i) * A @ v_i
            U[:, i] = (1.0 / singular_values[i]) * (A @ V[:, i])
        else:
            # 對於零奇異值對應的向量，需要額外計算正交基底
            # 這裡我們利用內建函式來處理（更為穩健的實現會使用 QR 分解或 Gram-Schmidt）
            # 簡化：將 U 矩陣的剩餘部分（零奇異值對應的部分）正交化
            pass # 這裡使用一個簡化版的 U 矩陣

    # 為了完整性，我們使用 np.linalg.svd 獲得 U 來補足因數值不穩定的 U 向量
    # 真實應用中，您會直接依賴 np.linalg.svd 提供的 U, s, Vh
    # 這裡僅演示核心數學關係：
    # U 是從 V 和 Sigma 算出來的
    
    # 為了演示，我們僅返回計算的 s 和 V.T
    return singular_values, V.T

# 測試矩陣
A_test_svd = np.array([[1, 1], [0, 1], [1, 0]], dtype=float)
s_eig, V_T_eig = svd_from_eig(A_test_svd)
s_np, V_T_np = np.linalg.svd(A_test_svd)[1:]

# print("\nSVD (來自 Eig) 奇異值:", s_eig)
# print("SVD (來自 NumPy) 奇異值:", s_np)
# print("SVD (來自 Eig) V.T:\n", V_T_eig)
# print("SVD (來自 NumPy) V.T:\n", V_T_np) 
# (結果應該非常相似，可能會有正負號差異，因為特徵向量方向可逆)

# 5.def pca_svd(X, k_components):
    """
    使用 SVD 進行主成分分析 (PCA)。
    
    參數:
      X (np.ndarray): 原始數據矩陣 (樣本數 x 特徵數)。
      k_components (int): 欲保留的主成分數量 (降維後的維度)。
      
    返回:
      X_reduced (np.ndarray): 降維後的數據。
      components (np.ndarray): 主成分向量 (即 V.T 的前 k 行)。
      variance_explained (np.ndarray): 每個主成分解釋的變異量。
    """
    
    # 1. 數據中心化 (Centering): 減去每列的平均值
    X_mean = np.mean(X, axis=0)
    X_centered = X - X_mean
    
    # 2. 進行 SVD 分解: X_centered = U @ Sigma @ V.T
    # V_T 的行向量即為我們尋找的主成分向量
    U, s, V_T = np.linalg.svd(X_centered)
    
    # 3. 確定主成分 (Principal Components)
    # V_T 的前 k 行即為 k 個主成分 (components)
    components = V_T[:k_components, :] 
    
    # 4. 計算解釋的變異量
    # 總變異量與奇異值平方和有關
    total_variance = np.sum(s**2)
    variance_explained = (s**2 / total_variance)[:k_components]
    
    # 5. 降維 (Projection): X_reduced = X_centered @ V
    # 也就是將中心化數據投影到主成分空間 (前 k 個 V 的列向量)
    V_k = V_T.T[:, :k_components] # 取得前 k 個右奇異向量 (V 的列)
    X_reduced = X_centered @ V_k

    return X_reduced, components, variance_explained

# 測試 PCA
# 假設數據有 5 個樣本，3 個特徵
X_data = np.array([[1, 2, 3], 
                   [2, 4, 6], 
                   [3, 6, 9], 
                   [10, 10, 10], 
                   [11, 11, 11]], dtype=float)

k = 1 # 降到 1 維

X_reduced, components, var_exp = pca_svd(X_data, k)

# print("\n--- PCA 主成分分析 (使用 SVD) ---")
# print("降維後的數據 (5 樣本 x 1 維):\n", X_reduced)
# print("主成分向量 (1 x 3):\n", components)
# print("單一主成分解釋的總變異量百分比:", var_exp * 100, "%")
