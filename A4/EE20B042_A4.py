'''
NAME:HARIHARAN P
ROLL NO.:EE20B042
EE2703 APL ASSIGNMENT 4: FOURIER APPROXIMATION
'''

# Importing the required modules	
from scipy import integrate
from pylab import *

#initialising the given two function
def e(x):
	return exp(x)

def coscos(x):
	return cos(cos(x))

#function to return the first n fourier coefficients of the function f using the inteegrate.quad 
def F_C(n,function):
    
    coeff= zeros(n)
    def fcos(x,k,f):
        return f(x)*cos(k*x)/pi
    def fsin(x,k,f):
        return f(x)*sin(k*x)/pi

    coeff[0] = integrate.quad(function,0,2*pi)[0]/(2*pi)
    
    for i in range(1,n):
        if(i%2==1):
            coeff[i] = integrate.quad(fcos,0,2*pi,
            args=(int(i/2)+1,function))[0]
        else:
            coeff[i] = integrate.quad(fsin,0,2*pi,
            args=(int(i/2),function))[0]
    return coeff 

#range of x values from -2π to 4π
x_tot = linspace(-2*pi,4*pi,1200)

#the corresponding function values for those x values
exp_val = e(x_tot)
coscos_val = coscos(x_tot)

#making exp() and coscos(x) periodic 
p_exp_val = []
p_exp_val[0:400] = exp_val[400:800]
p_exp_val[400:800] = exp_val[400:800]
p_exp_val[800:1200] = exp_val[400:800]

p_coscos_val = []
p_coscos_val[0:400] = coscos_val[400:800]
p_coscos_val[400:800] = coscos_val[400:800]
p_coscos_val[800:1200] = coscos_val[400:800]

#getting the functions' coeff using integration method
e_coeff=F_C(51,e)
coscos_coeff=F_C(51,coscos)

x=linspace(0,2*pi,401)
x=x[:-1] # drop last term to have a proper periodic integral

#getting the actual vlues of the function
e_b=e(x)
coscos_b=coscos(x)

A=zeros((400,51)) # initialising A
A[:,0]=1 # col 1 is all ones
for k in range(1,26):
    A[:,2*k-1]=cos(k*x) # cos(kx) column
    A[:,2*k]=sin(k*x) # sin(kx) column

e_c=lstsq(A,e_b,rcond=None)[0] 
coscos_c=lstsq(A,coscos_b,rcond=None)[0] 
# The difference between the actual and predicted values are found
e_diff= abs(e_coeff - e_c)
coscos_diff = abs(coscos_coeff - coscos_c)
# The deviation between the actual and predicted values are found
e_dev = e_diff.max()
coscos_dev = coscos_diff.max()
# The estimated values of the function using matrix method 
estimated_e=dot(A,e_c)
estimated_coscos=dot(A,coscos_c)

print("The deviation between coefficients for exp() is ",e_dev)
print("The deviation between coefficients for coscos() is ",coscos_dev)

#plot all the figures asked
figure(1)
semilogy(x_tot,exp_val,label='True')
semilogy(x_tot,p_exp_val,label='Periodic Extension')
semilogy(x,estimated_e,'og',label='Predicted')
title("Figure 1: exp(x) function")
xlabel(r'$x\rightarrow$')
ylabel(r'$e^x\rightarrow$')
grid(True)
legend()

figure(2,figsize=(10, 8), dpi=80)
plot(x_tot,coscos_val,label='True')
semilogy(x_tot,p_coscos_val,label='Periodic Extension')
semilogy(x,estimated_coscos,'og',label='Predicted')
title("Figure 2: cos(cos(x)) function")
xlabel(r'$x\rightarrow$')
ylabel(r'$cos(cos(x))\rightarrow$')
grid(True)
legend()

n=list(range(51))

figure(3)
semilogy(n,abs(e_coeff),'or',label='True')
semilogy(n,abs(e_c),'og',label='Predicted')
title("Figure 3: coeff of exp(x) in semilog")
xlabel(r'$n\rightarrow$')
ylabel(r'$coefficients of e^x\rightarrow$')
grid(True)
legend()

figure(4)
loglog(n,abs(e_coeff),'or',label='True')
loglog(n,abs(e_c),'og',label='Predicted')
title("Figure 4: coeff of exp(x) in loglog")
xlabel(r'$n\rightarrow$')
ylabel(r'$coefficients of e^x\rightarrow$')
grid(True)
legend()

figure(5)
semilogy(n,abs(coscos_coeff),'or',label='True')
semilogy(n,abs(coscos_c),'og',label='Predicted')
title("Figure 5: coeff of cos(cos(x)) in semilog")
xlabel(r'$n\rightarrow$')
ylabel(r'$coefficients of coscos(X)\rightarrow$')
grid(True)
legend()

figure(6)
loglog(n,abs(coscos_coeff),'or',label='True')
loglog(n,abs(coscos_c),'og',label='Predicted')
title("Figure 6: coeff of cos(cos(x)) in loglog")
xlabel(r'$n\rightarrow$',size=15)
ylabel(r'$coefficients of coscos(x)\rightarrow$')
grid(True)
legend()

show()