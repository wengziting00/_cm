#1. 定義基本幾何物件
import math

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y
        
    def __repr__(self):
        return f"Point({self.x}, {self.y})"


class Line:
    # ax + by + c = 0
    def __init__(self, a, b, c):
        self.a = a
        self.b = b
        self.c = c

    @classmethod
    def from_points(cls, p1, p2):
        a = p2.y - p1.y
        b = p1.x - p2.x
        c = -(a * p1.x + b * p1.y)
        return cls(a, b, c)
    
    def __repr__(self):
        return f"{self.a}x + {self.b}y + {self.c} = 0"


class Circle:
    def __init__(self, center: Point, r):
        self.center = center
        self.r = r
        
    def __repr__(self):
        return f"Circle(center={self.center}, r={self.r})"
    
#2.求兩直線交點
def line_intersection(L1, L2):
    D = L1.a * L2.b - L2.a * L1.b
    if abs(D) < 1e-9:
        return None  # 平行或重合
    
    x = (L1.b * L2.c - L2.b * L1.c) / D
    y = (L2.a * L1.c - L1.a * L2.c) / D
    return Point(x, y)

#3.直線與圓交點
def line_circle_intersection(line, circle):
    a, b, c = line.a, line.b, line.c
    x0, y0, r = circle.center.x, circle.center.y, circle.r

    # 平移圓心到原點（簡化）
    c2 = a*(x0) + b*(y0) + c

    # 解: (ax + by + c2)^2 = r^2(a^2 + b^2)
    A = a*a + b*b
    B = 2*a*c2
    C = c2*c2 - r*r*A

    disc = B*B - 4*A*C
    if disc < 0:
        return []   # 無交點
    elif abs(disc) < 1e-9:
        t = -B / (2*A)
        x = x0 - a*t
        y = y0 - b*t
        return [Point(x,y)]
    else:
        sqrtD = math.sqrt(disc)
        t1 = (-B + sqrtD) / (2*A)
        t2 = (-B - sqrtD) / (2*A)
        return [
            Point(x0 - a*t1, y0 - b*t1),
            Point(x0 - a*t2, y0 - b*t2)
        ]
 #4. 兩圓交點
def circle_intersection(c1, c2):
    x1, y1, r1 = c1.center.x, c1.center.y, c1.r
    x2, y2, r2 = c2.center.x, c2.center.y, c2.r

    dx = x2 - x1
    dy = y2 - y1
    d = math.hypot(dx, dy)

    if d > r1 + r2 or d < abs(r1 - r2) or d == 0:
        return []  # 無交點或圓重合
    
    a = (r1*r1 - r2*r2 + d*d) / (2*d)
    h = math.sqrt(r1*r1 - a*a)

    xm = x1 + a * dx / d
    ym = y1 + a * dy / d

    rx = -dy * (h/d)
    ry = dx * (h/d)

    return [
        Point(xm + rx, ym + ry),
        Point(xm - rx, ym - ry)
    ]
#5. 求線外一點的垂足（投影）
def foot_of_perpendicular(line, p):
    a, b, c = line.a, line.b, line.c
    denom = a*a + b*b
    d = (a*p.x + b*p.y + c)

    x = p.x - a * d / denom
    y = p.y - b * d / denom
    return Point(x, y)
#6.利用垂足驗證畢氏定理
def dist2(p, q):
    return (p.x - q.x)**2 + (p.y - q.y)**2

def verify_pythagoras(A, B, C):
    return abs(dist2(A,B) - (dist2(A,C) + dist2(C,B))) < 1e-6
#7. 定義三角形物件
class Triangle:
    def __init__(self, A, B, C):
        self.A = A
        self.B = B
        self.C = C
    
    def vertices(self):
        return (self.A, self.B, self.C)
#8.幾何物件的平移、縮放、旋轉
#  平移
def translate_point(p, dx, dy):
    return Point(p.x + dx, p.y + dy)

#縮放
def scale_point(p, sx, sy):
    return Point(p.x * sx, p.y * sy)

#旋轉
def rotate_point(p, theta):
    c = math.cos(theta)
    s = math.sin(theta)
    return Point(p.x*c - p.y*s, p.x*s + p.y*c)
