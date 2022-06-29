'''
Name: Hariharan P
Roll No.:EE20B0442
EE2703 APL ENDSEM: HALF WAVE DIPOLE ANTENNAS
'''
# Import the required libraries.
from pylab import *

#Initialise the required variables
n=4  #no. of sections 
dz=0.5/n # length of section
a=0.01 #radius
k=pi # wavenumber
mu=(4*pi*(10**-7))/(4*pi) #mu0/4pi
muo=4*pi*(10**-7) #mu0

#qn1: Initialising I and J current vectors and z and u vectors
# z is initiaised as per z = i×dz, −N ≤ i ≤ N and u without first,middle,last element of z 
I=zeros(2*n+1) #current vector using equation given
J=zeros(2*n-2) #current vector to be computed
z=array(arange(-0.5,0.5+dz,dz)) #vector z 
u=array(append(arange(-0.5+dz,0,dz),arange(dz,0.5,dz))) #vector u

#qn2: a function to return m matrix
def getM():
  m=identity(2*n-2)
  m=m/(2*pi*a) #m matrix is 1/(2*pi*a) times indentity matrix
  return m

#qn3:calculation of P and Pb matrices as per the formulae given
(zj,zi)=meshgrid(z,z)
(uj,ui)=meshgrid(u,u)
Rz=array((a**2+(zi-zj)**2)**0.5) #distance from source to observer for all current elements
Ru=array((a**2+(ui-uj)**2)**0.5) #distance from source to observer for unknown current elements

P=zeros((2*n-2,2*n-2),dtype=complex_) #initialising P and Pb matrices as zeros
Pb=zeros(2*n-2,dtype=complex_)

RiN=delete(Rz[:,n],[0,n,2*n]) #slicing the middle row from Rz and removing first,middle,last elements to make RiN vector
def epower(x):  # function to return complex value if the input to e is complex(eular's form)
  return (complex(cos(x),sin(x)))

epower=vectorize(epower) #vectorising the function
epower1 =epower(-1*Ru*k)  
epower2 =epower(-1*k*RiN)
epower3 =epower(pi/2)
P=(mu*dz*epower1)/(Ru)   #preparing P
Pb=(mu*dz*epower2)/(RiN)  #preparing Pb

#q4: computing Q and Qb matrices
Q=zeros((2*n-2,2*n-2),dtype=complex_) #initialising Q and Qb
Qb=zeros(2*n-2,dtype=complex_)
Q=P*(a/muo)*(((k*epower3)/Ru)+(1/(Ru**2)))  #preparing Q
Qb=Pb*(a/muo)*(((k*epower3)/RiN)+(1/(RiN**2))) #preparing Qb

#qn5: computing J current vector using linalg inverse on (m-Q) matrix and multiplying it with Qb*Im
m=getM() #calling getM function to get m matrix
J=dot(linalg.inv(m-Q),(Qb*1.0)) #computing J matrix which is for unknown current
J=insert(J,0,0)  #include known currents
J=append(J,0)
J=insert(J,n,1)
z1=z[:n]
z2=z[n:2*n+1]
I1=array(sin(k*(0.5+z1))) #current for negative values of z using the given direct formula
I2=array(sin(k*(0.5-z2))) #current for positive values of z using the given direct formula
I=concatenate((I1,I2)) #concatenating to form a single vector of known currents

#plotting computed current and current found using direct formula VS z
plot(z,abs(J),label='computed current',color='r')
plot(z,I,label='equation current',color='g')
title("current VS distance")
xlabel(r'$distance\rightarrow$')
ylabel(r'$current\rightarrow$')
legend()
grid(True)
show()
'''
#Printing all matrices used
print("z:\n",(z).round(2))
print("u:\n",(u).round(2))
print("I:\n",(I).round(2))
print("abs(J * 10**3) :\n",abs(J*((10)**3)).round(2))
print("m:\n",(m).round(2))
print("Rz:\n",(Rz).round(2))
print("Ru:\n",(Ru).round(2))
print("P*((10)**8):\n",(P*((10)**8)).round(2))
print("Pb*((10)**8): \n",(Pb*((10)**8)).round(2))
print("Q*((10)**4):\n",(Q*((10)**4)).round(2))
print("Qb*((10)**4):\n",(Qb*((10)**4)).round(2))
'''