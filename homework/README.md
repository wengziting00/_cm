## 習題 1 : 請用程式驗證微積分基本定理

[習題](https://github.com/wengziting00/_cm/blob/main/homework/hw1/hw1.py)
[說明](https://github.com/wengziting00/_cm/blob/main/homework/hw1/README.md)
設定一個非常小的數 ℎ=0.00001 h=0.00001，用來模擬微分與積分的極小變化量。接著利用「差分法」定義 df(f, x)，以 (f(x+h)−f(x))/h 來近似函數在 x 處的導數；再用「黎曼和」的概念定義 integral(f, a, b)， 從 a 開始以步長 h 累加每一小段的面積 f(x)⋅h，近似計算定積分。之後在 theorem1(f, x) 中，先把函數 f 從 0 積分到 x，再對這個積分結果做微分，理論上應該回到原函數f(x)， 這正是在數值方式下驗證微積分基本定理。程式最後以 𝑓(𝑥)=x^3為例，計算在 x=2 處的導數近似值與從 0 到 2 的定積分近似值，並檢查「先積分再微分」的結果是否與原函數在該點的值足夠接近（誤差小於 0.01）
