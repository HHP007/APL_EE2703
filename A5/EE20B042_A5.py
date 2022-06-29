'''
Name: Hariharan P
Roll No.:EE20B0442
EE2703 APL ASSIGNMENT 5: THE RESISTOR PROBLEM
'''

# Import the required libraries.	
from pylab import *
import sys
import mpl_toolkits.mplot3d.axes3d as p3

#getting inputs - default/user
if(len(sys.argv)==5):
    Nx=int(sys.argv[1])
    Ny=int(sys.argv[2])
    radius=int(sys.argv[3])  
    Niter=int(sys.argv[4])
    print("user provided parameters are being used")
else:
    Nx=25
    Ny=25 
    radius=8 
    Niter=1500 
    print("Using default values")

#Allocate the potential array and initialize it
x,y=linspace(-0.5,0.5,num=Nx,dtype=float),linspace(-0.5,0.5,num=Ny,dtype=float)
X,Y = meshgrid(x,-y)
ii=where(X**2+Y**2<(0.35)**2)
phi=zeros((Nx,Ny),dtype = float)
phi[ii]=1.0

figure(1)
plot((ii[0]-Nx/2)/Nx,(ii[1]-Ny/2)/Ny,'ro',label="V = 1")
title("Countor Plot of Initial Potential ")
xlim(-0.5,0.5)
ylim(-0.5,0.5)
xlabel(r'$X\rightarrow$')
ylabel(r'$Y\rightarrow$')
grid(True)
legend()
show()
#performing the iteration and calculating the error in the potential
errors = empty(Niter)
for k in range(Niter):
	oldphi = phi.copy()
	phi[1:-1,1:-1] = 0.25*(phi[1:-1,0:-2] + phi[1:-1,2:]
     + phi[0:-2,1:-1] + phi[2:,1:-1])
	phi[1:-1,0] = phi[1:-1,1]
	phi[1:-1,-1] = phi[1:-1,-2]
	phi[0,1:-1] = phi[1,1:-1]
	phi[ii] = 1.0
	errors[k]=(abs(phi-oldphi)).max();
#getting the best fit
def best_fit(y,Niter,lastn=0):
    log_err = log(y)[-lastn:]
    X = vstack([(np.arange(Niter)+1)[-lastn:],ones(log_err.shape)]).T
    log_err = reshape(log_err,(1,log_err.shape[0])).T
    return lstsq(X, log_err,rcond=None)[0]
#computing the Cumulative error
def cummu_error(a,b,Niter):
    return -a/b*exp(b*(Niter+0.5))

b,a = best_fit(errors,Niter)
b_,a_ = best_fit(errors,Niter,500)

figure(2)
title("Best fit for error on a loglog scale")
xlabel("No of iterations")
ylabel("Error")
xnew = asarray(range(Niter))+1
loglog(xnew,errors)
loglog(xnew[::100],exp(a+b*asarray(range(Niter)))[::100],'ro')
loglog(xnew[::100],exp(a_+b_*asarray(range(Niter)))[::100],'go')
legend(["errors","fit1","fit2"])
show()

iteration=arange(100,Niter+1,100)
figure(3)
grid(True)
title(r'Plot of Cumulative Error values On a loglog scale')
loglog(iteration,abs(cummu_error(a_,b_,iteration)),'ro')
xlabel("iterations")
ylabel("Net  maximum error")
show()
#exponent part of the error values
c_approx = lstsq(c_[ones(Niter),arange(Niter)]
,log(errors),rcond=None)
a, b = exp(c_approx[0][0]), c_approx[0][1]
print(" A and B : ",a,b)
c_500 = lstsq(c_[ones(Niter-500),arange(500,Niter)]
,log(errors[500:]),rcond=None)
a__,b__ = exp(c_500[0][0]),c_500[0][1]
print("A and B for the iterations after 500 : ",a__,b__)

# The current densities are calculated 
Jx = zeros((Ny, Nx))
Jy = zeros((Ny, Nx))
Jx[1:-1, 1:-1] = 0.5*(phi[1:-1, 0:-2] - phi[1:-1, 2:])
Jy[1:-1, 1:-1] = 0.5*(phi[2:, 1:-1] - phi[0:-2, 1:-1])

n = arange(Niter)
#error vs iteration in semilog
figure(4)
semilogy(n,errors)
title("Error versus iteration")
xlabel(r'$Iteration\rightarrow$',size=15)
ylabel(r'$Error\rightarrow$',size=15)
grid(True)

#error vs iteration above 500 in semilog 
figure(5)
semilogy(n[500:],errors[500:])
title("Error versus iteration above 500")
xlabel(r'$Iteration\rightarrow$',size=15)
ylabel(r'$Error\rightarrow$',size=15)
grid(True)

#actual and expected error (above 500 iterations) in semilog
figure(6)
semilogy(n[500:],errors[500:],label="Actual")
semilogy(n[500:],a__*exp(b__*n[500:]),label="Expected")
title("Expected versus actual error (>500 iterations)")
xlabel(r'$Iteration\rightarrow$',size=15)
ylabel(r'$Error\rightarrow$',size=15)
grid(True)
legend()

#the actual and expected error in semilog
figure(7)
semilogy(n,errors,label="Actual")
semilogy(n,a*exp(b*n),label="Expected")
title("Expected versus actual error ")
xlabel(r'$Iteration\rightarrow$',size=15)
ylabel(r'$Error\rightarrow$',size=15)
grid(True)
legend()

#the contour of phi 
figure(8)
contourf(X,Y,phi)
plot((ii[0]-Nx/2)/Nx,(ii[1]-Ny/2)/Ny,'ro',label="V = 1")
title("Contour plot of potential")
xlabel(r'$X\rightarrow$')
ylabel(r'$Y\rightarrow$')
colorbar()
grid(True)
legend()

#the surface plots of phi 
fig1=figure(9)
ax=p3.Axes3D(fig1) 
title("The 3-D surface plot of the potential")
xlabel(r'$X\rightarrow$')
ylabel(r'$Y\rightarrow$')
surf = ax.plot_surface(X, Y, phi, rstride=1, cstride=1, cmap=cm.jet)
fig1.colorbar(surf)

#the vector plot of current flow
figure(10)
quiver(X,Y,Jx,Jy)
plot((ii[0]-Nx/2)/Nx,(ii[1]-Ny/2)/Ny,'ro')
title("The vector plot of the current flow")
xlabel(r'$X\rightarrow$')
ylabel(r'$Y\rightarrow$')
show()
