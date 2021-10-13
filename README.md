import numpy as np
import matplotlib as plt
import math
import scipy
from scipy.integrate import odeint

k = (2800000)
w = (100000)
zta = (.05)
f0 = (100000)
wn = math.sqrt((k*32.2)/w)
dt = (0.005)
wd = ((math.sqrt(1 - zta**2))*wn)
bt = (zta*wn)
f1 = (wd*dt)
f2 = (bt*dt)
f3 = (k*wd*dt)
f4 = (math.exp(-f2))
f5 = (math.sin(f1))
f6 = (math.cos(f1))
f7 = (math.cos(-f1))

a = (((f4*(wd**2 - bt**2)/wn**2 - f2)*f5 -((2*wd*bt)/wn**2 + f1)*f6 + ((2*wd*bt)/wn**2))/f3) 
b = (((f4*(wd**2 - bt**2)/wn**2 - f2)*f5 +((2*wd*bt)/wn**2)*f6 - ((2*wd*bt)/wn**2))/f3)
c = (f4*f6 + (bt/wd)*f5)
d = ((f4*f5)/wd)
a1 = ((f4*((bt + wn**2*dt)*f5 + wd*f6) - wd)/f3)
b1 = (((-f4)*(bt*f5 + wd*f6) + wd)/f3)
c1 = ((wn**2/wd)*(-f4)*f5)
d1 = (f4*f6-(bt/wd)*f5)
print(a)
print(b)
print(c)
print(d) 
print(a1)
print(b1)
print(c1)
print(d1)

#initial conditions
x_0 = [f0, f0, 0, 0]

dti = np.arange(0, 0.5, dt)

def tank(x, dti):
    dx1_dt = f0*(1 -(dti - dt)/0.05)
    dx2_dt = f0*(1 - dti/0.05)
    dx3_dt = a*x1 + b*x2 + c*x3 + d*x4
    dx4_dt = a1*x1 + b1*x2 + c1*x3 + d1*x4
    dx_dt = [dx1_dt, dx2_dt, dx3_dt, dx4_dt]
    return dx_dt

x_results = odeint(tank, x_0, dti)
p = []
p1 = []
d = []
v = []

for i in range(0, len(x_results)):
    #displacement
     d_value = x_results[i, 0]
     d_value.append(d)
    

    #velocity
     v_value = x_results[i, 0]
     v_value.append(v)
    
