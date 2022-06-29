'''
Name: Hariharan P
Roll no.: EE20B042
APL EE2703 A7 ASSIGNMENT 
'''

from sympy import *
import scipy.signal as sp
import pylab as plt
s = symbols('s')

def lowpass(R1,R2,C1,C2,G,Vi):

	A = Matrix([[0,0,1,-1/G],
		[-1/(1+s*R2*C2),1,0,0],
		[0,-G,G,1],
		[-1/R1-1/R2-s*C1,1/R2,0,s*C1]])
	b = Matrix([0,0,0,-Vi/R1])
	V = A.inv()*b
	return A,b,V


def highpass(R1,R3,C1,C2,G,Vi):

    A = Matrix([[0,-1,0,1/G],
        [s*C2*R3/(s*C2*R3+1),0,-1,0],
        [0,G,-G,1],
        [-s*C2-1/R1-s*C1,0,s*C2,1/R1]])
    b = Matrix([0,0,0,-Vi*s*C1])
    V = A.inv()*b
    return A,b,V

 	 
def sympyToTrFn(Y):
	
    Y = expand(simplify(Y))
    
    num,den = fraction(Y)
    num, den = [[float(i) for i in Poly(j, s).all_coeffs()] for j in Y.as_numer_denom()]
  
    H = sp.lti(num,den)
    return H



A,b,V=lowpass(10000,10000,1e-9,1e-9,1.586,1)
Vout=V[3]
H = sympyToTrFn(Vout)
ww=plt.logspace(0,8,801)
ss=1j*ww
hf=lambdify(s,Vout,'numpy')
v=hf(ss)

plt.title("Low pass Magnitude response")
plt.xlabel(r'$\omega\rightarrow$')
plt.ylabel(r'$|H(j\omega)|\rightarrow$')
plt.loglog(ww,abs(v),lw=2)
plt.grid(True)
plt.show() 
  
As1,Bs1,Vs1=lowpass(10000,10000,1e-9,1e-9,1.586,1/s)
Vout1=Vs1[3]
H1 = sympyToTrFn(Vout1)
t,y1 = sp.impulse(H1,None,plt.linspace(0,5e-3,10000))
plt.plot(t,y1)
plt.title("Step Response for low pass filter")
plt.xlabel(r'$t\rightarrow$')
plt.ylabel(r'$V_o(t)\rightarrow$')
plt.grid(True)
plt.show() 


A,b,V2=highpass(10000,10000,1e-9,1e-9,1.586,1)
Vout2=V2[3]
H2 = sympyToTrFn(Vout2)
ww=plt.logspace(0,8,801)
ss=1j*ww
hf=lambdify(s,Vout2,"numpy")
v=hf(ss)

plt.title("High pass Magnitude response")
plt.xlabel(r'$\omega\rightarrow$')
plt.ylabel(r'$|H(j\omega)|\rightarrow$')
plt.loglog(ww,abs(v),lw=2)
plt.grid(True)
plt.show() 
  
As2,Bs2,Vs2=highpass(10000,10000,1e-9,1e-9,1.586,1/s)
Vout4=Vs2[3]
H3 = sympyToTrFn(Vout2)
t,y2 = sp.impulse(H3,None,plt.linspace(0,5e-3,10000))
plt.plot(t,y2)
plt.title("Step Response for high pass filter")
plt.xlabel(r'$t\rightarrow$')
plt.ylabel(r'$V_o(t)\rightarrow$')
plt.grid(True)
plt.show() 

t=plt.linspace(0,1e-2,100000)
def inp1(t):
    return (plt.sin(2000*plt.pi*t)+plt.cos(2e6*plt.pi*t))
def inp2(t):
    return plt.sin(1e7*t)*plt.exp(-1000*t) 
def inp3(t):
    return plt.sin(1e3*t)*plt.exp(-1000*t) 
def inp_response(Y,title,inp,tmax=5e-3):
    H= sympyToTrFn(Y)
    t = plt.linspace(0,tmax,100000)
    t,y,_ = sp.lsim(H,inp(t),t)
    plt.title(title)
    plt.xlabel("t")
    plt.ylabel("Vo")
    plt.plot(t,y)
    plt.show()
    return


inp_response(Vout,"Response of Low pass Filter to sum of sinusoids"
,inp1)
inp_response(Vout2,"Response of High pass Filter to sum of sinusoids"
,inp1,1e-2)

inp_response(Vout2,
"Response of high pass Filter to Damped High Frequency sinusoid"
,inp2,1e-2)
inp_response(Vout2,
"Response of high pass Filter to Damped low Frequency sinusoid"
,inp3,1e-2)