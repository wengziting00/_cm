#1

import numpy as np

def dft(x):

    N = len(x)
    n = np.arange(N) 
    k = n.reshape((N, 1)) 

    M = np.exp(-2j * np.pi * k * n / N)
     
    X = M.dot(x)
    
    return X

x_test = np.array([1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0])
X_dft = dft(x_test)

#2
def idft(X):

    N = len(X)
    k = np.arange(N) 
    n = k.reshape((N, 1)) 
    M = np.exp(2j * np.pi * k * n / N)
    
    x = M.dot(X)
    
    x = x / N
    
    return x

x_idft = idft(X_dft)


#3
def verify_dft_inverse(f):

    
    print("\n--- 驗證 DFT/IDFT 可逆性 ---")
    
    F = dft(f)
    print(f"原始信號長度 N={len(f)}")
   
    
   
    f_recovered = idft(F)
    f_recovered_real = f_recovered.real
    
    error = np.linalg.norm(f - f_recovered_real)
    
    print(f"原始信號 f:\n {f}")
    print(f"逆變換 f' 的實部:\n {f_recovered_real.round(decimals=10)}")
    print(f"\n原始信號與逆變換結果的 L2 誤差: {error:.2e}")
    
    if error < 1e-9:
        print("結論: 誤差極小，驗證成功，DFT 和 IDFT 互為逆運算。")
    else:
        print("結論: 誤差較大，驗證失敗。")


N_samples = 64
t = np.linspace(0, 1, N_samples, endpoint=False)
f_test = 2 * np.sin(2 * np.pi * 5 * t) + 1.5 * np.cos(2 * np.pi * 12 * t) 

verify_dft_inverse(f_test)

