import math as m
import numpy as np
pie=m.pi
c=3*10**8
Er=eps=2.2
fr=12*10**9
h_mm=1.575
h_m=h_mm*0.001
w=15
μ0 =4*pie*1e-7

h=h_mm*0.1
epsiloneff=((eps+1)/2)+((eps-1)/2)*(1+10*(h/w))**(-0.5)
F=8.791*10**9/(fr*m.sqrt(eps))
a=F/m.sqrt(1+((2*h)/(pie*eps*F))*(1.7726+m.log2(pie*F*0.5/h)))
print("Patch Radius in cm : ",a)

h_cm = h_mm/10.0
    
F = 8.791e9/(fr * np.sqrt(Er))

a_cm = np.power((1+ ((2*h_cm/(np.pi*Er*F))* (np.log(0.5*np.pi*F/h_cm) + 1.7726))), -0.5)
a_cm = a_cm * F
a = a_cm * 10.0

ae_cm = np.log(np.pi*a_cm*0.5/h_cm) + 1.7726
ae_cm = 2.0*h_cm*ae_cm/(np.pi*a_cm*Er)
ae_cm = np.sqrt(1.0 + ae_cm)
ae_cm = a_cm * ae_cm

ae = ae_cm * 10.0
print('Effective Radius of Circular Patch =',ae,'mm') 

lamda=c/fr
print("ae/lamda for Grad",ae*0.001/lamda)
Grad=1e-3
Gd=2.39*0.0019/(μ0*h_m*fr )
print("G due to dielectric",Gd)
Gc=2.39*pie*m.pow(pie*fr*μ0, -1.5)/(4*h_m*h_m*m.sqrt(58.14e6))
print("G due to Conduction loss",Gc)
Gt=Grad+Gd+Gc
print("Total conductance",Gt)
Z0=50
Zin=1/Gt
Zl=m.sqrt(Z0*Zin)
print("Load Impedance(Patch)",Zl)
Z=Zl
A=((Z/60)*(m.sqrt(0.5*eps+0.5)))+((eps-1)/(eps+1))*(0.23+(0.11/eps))
x=8*m.exp(A)/(m.exp(2*A)-2)
#print("w1/h : ",x)
print("Width of quarter wave transformer w1",x*h_mm,'mm')
Z=50
A=((Z/60)*(m.sqrt(0.5*eps+0.5)))+((eps-1)/(eps+1))*(0.23+(0.11/eps))
x=8*m.exp(A)/(m.exp(2*A)-2)
#print("w1/h : ",x)
print("Width of feed line w1",x*h_mm,'mm')

print(lamda)


B=377*m.pi/(2*Z*m.sqrt(eps))
y=(2/m.pi)*(B-1-m.log(2*B-1)+(((eps-1)/(eps+1))*(m.log(B-1)+0.39-(0.61/eps))))
print("w1/h : ",y)
print("Width of feed line w1",y*h_mm,'mm')



lamda=c/fr
w=lamda/(2*m.pow((eps+1)/2, -0.5))
print("Optimum width of feed line",w)

h=1.575
z0=50
w=(377*h)/(z0*m.sqrt(eps))
print("Width of feedline ",w)
