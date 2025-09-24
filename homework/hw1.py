h = 0.00001 # 用非常小的質來模擬微分與積分的計算

# 使用差分法來模擬微分
def df(f, x): # 定義df 參數f x 
    return (f(x+h)-f(x))/h 
# f(x + h)：函數在x+h的值
#f(x)：函數在x處的值
#h = 0.0001


def integral(f, a, b): #定義函數integral 參數fab
    x = a #積分起點a
    area = 0 #起始0
    while x < b: #x還沒到b時
        area += f(x)*h #f(x)*h累加進area
        x+=h #每次增加
    return area #回傳累加後的總面積

def theorem1(f, x): #定義
    r = df(lambda x:integral(f, 0, x), x) 
    #lambda x:integral(f, 0, x)唯一個函數傳給df() 對x處做微分計算
    print('r=', r, 'f(x)=', f(x)) #顯示微分和原始函數值
    print('abs(r-f(x))<0.01 = ', abs(r-f(x))<0.01) #顯示兩者是否接近
    assert abs(r-f(x))<0.01 #大於0.01失敗

def f(x): #f(x)=x3
    return x**3

print('df(f, 2)=', df(f, 2)) #輸出f(x)在x=2的導數近似值
print('integral(f, 0, 2)=', integral(f, 0, 2)) #輸出f(x)從x=0到x=2的定積分近似值

theorem1(f, 2) #驗證
