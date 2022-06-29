'''
Name: Hariharan P
Roll No.: EE20B042
APL ASSIGNMENT 6 :THE LAPLACE TRANSFORM
'''
# Import the necessary libraries
from pylab import *
import scipy.signal as sp
#Q1 and Q2
def spring_tf(f,decay):
    p= polymul([1,-2*decay,f*f + decay*decay],[1.0,0,2.25])
    return sp.lti([1,-1*decay],p)

t1,x1 = sp.impulse(spring_tf(1.5,-0.5),None,linspace(0,50,501))
t2,x2 = sp.impulse(spring_tf(1.5,-0.05),None,linspace(0,50,501))

plot(t1,x1)
title("The solution x(t) for delay=-0.5")
xlabel(r'$t\rightarrow$')
ylabel(r'$x(t)\rightarrow$')
grid(True)
show()

plot(t2,x2)
title("The solution x(t) for delay=-0.05")
xlabel(r'$t\rightarrow$')
ylabel(r'$x(t)\rightarrow$')
grid(True)
show()
#Q3
title("Forced Damping Oscillator with Different frequencies")
xlabel(r'$t\rightarrow$')
ylabel(r'$x(t)\rightarrow$')
details = []
for w in arange(1.4,1.6,0.05):
    H = sp.lti([1],[1,0,2.25])
    t = linspace(0,50,501)
    func = cos(w*t)*exp(-0.05*t)*(t>0)
    t,y,svec = sp.lsim(H,func,t)
    plot(t,y)
    details.append("frequency:" + str(w))
grid(True)
legend(details)
show()
#Q4
X3 = sp.lti([1,0,2],[1,0,3,0])
Y3 = sp.lti([2],[1,0,3,0])	
t3,x3 = sp.impulse(X3,None,linspace(0,20,501))
t3,y3 = sp.impulse(Y3,None,linspace(0,20,501))
plot(t3,x3,label='x(t)')
plot(t3,y3,label='y(t)')
title("Coupled Oscillations x(t),y(t)")
xlabel(r'$t\rightarrow$')
ylabel(r'$functions\rightarrow$')
legend(loc = 'upper right')
grid(True)
show()
#Q5
denom = poly1d([1e-12,1e-4,1])
H = sp.lti([1],denom)
w,S,phi = H.bode()
semilogx(w,S)
title("Magnitude Bode plot")
xlabel(r'$\omega\rightarrow$')
ylabel(r'$20\log|H(j\omega)|\rightarrow$')
grid(True)
show()

semilogx(w,phi)
title("Phase Bode plot")
xlabel(r'$\omega\rightarrow$')
ylabel(r'$\angle H(j\omega)\rightarrow$')
grid(True)
show()
#Q6
vi_l = cos(1e3*arange(0,30e-3,1e-7)) - cos(1e6*arange(0,30e-3,1e-7))
vi_s = cos(1e3*arange(0,30e-6,1e-7)) - cos(1e6*arange(0,30e-6,1e-7))
t4,vo_s,svec1 = sp.lsim(H,vi_s,arange(0,30e-6,1e-7))
t5,vo_l,svec2 = sp.lsim(H,vi_l,arange(0,30e-3,1e-7))

plot(t4,vo_s)
title("The Output Voltage for small time interval")
xlabel(r'$t\rightarrow$')
ylabel(r'$V_o(t)\rightarrow$')
grid(True)
show()

plot(t5,vo_l)
title("The Output Voltage for large time interval")
xlabel(r'$t\rightarrow$')
ylabel(r'$V_o(t)\rightarrow$')
grid(True)
show()

