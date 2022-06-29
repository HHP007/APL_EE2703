'''
NAME: HARIHARAN P
ROLL NO.: EE20B042
EE2703 ASSIGNMENT3
'''

#importing necessay libraries

from pylab import *
import scipy.special as special

givendata=loadtxt("fitting.dat")# loading the data file generated using given generte data py file

time=givendata[:,0] # extracting the first column which is the time

def g(t,A,B):  # defining the given function
    return A*special.jn(2,t)+B*t

act_f=g(time,1.05,-0.105) #the actual values
stdev=logspace(-1,-3,9)

x,y=givendata.shape #x has no. of datas
# printing nine other noisy plots + the actual plot.
figure(0)
title("Q4:Data to be fitted to theory",size=15)	
xlabel('time',size=15)
ylabel('g()+noise',size=15)
plot(time,act_f,label="actual Value",color='black',linewidth=4)
for i in range(1,10):
	plot(time,givendata[:,i],label="stdev=%.4f"%stdev[i-1])
grid(True)
legend()

# printing the error plot of the first column of data vs the actual plot.
col2=givendata[:,1] #extracting the 2nd column that has the first set of data
figure(1)
title("Q5:Data with Error Bars",size=15)	
xlabel('time',size=15)
ylabel('data points',size=15)
plot(time,act_f,label="f(t)",color='black',linewidth=4)
errorbar(time[::5],col2[::5],stdev[0],fmt='ro',label='errorbar/noise')
grid(True)
legend()

#generating M matrix mentioned from two column vectors
col1=zeros((x,1))
for j in range(x):
	col1[j] = (special.jn(2,time[j]))
M=c_[col1,time]
p = array([1.05,-0.105])
G = dot(M,p)
#part of Q6
if array_equal(G,g(time,1.05,-0.105)):
	print("generated matrix matches with the actual values")
else:
	print("generated matrix does not match with the actual values")	

# calculating the mean squared error for the the set of values given for A and B

Ai = [0.1*i for i in range(21)]
Bi = [0.01*i-0.2 for i in range(21)]
Err = zeros((21,21))
for i in range(21):
	for j in range(21):
		for k in range(x):
			Err[i][j] += ((col2[k]-g(time[k],Ai[i],Bi[j]))**2)
Err=Err/x


#  finding the best estimate of A and B and calcualting its error
Err_A = zeros((9,1))
Err_B = zeros((9,1))

for i in range(1,10):
	best_est_AB = linalg.lstsq(M,givendata[:,i],rcond=None)
	Err_A[i-1] = abs(best_est_AB[0][0]-1.05)
	Err_B[i-1] = abs(best_est_AB[0][1]+0.105)

#X AND Y COORD AS MESHGRID 
X,Y = meshgrid(Ai,Bi)
# printing the contour of mean squared error for different values of A and B.
figure(2)
Contour=contour(X,Y,Err,15)
title("Q8:Contour Plot of Err_ij",size=15)
plot(1.05,-0.105,'b*')
annotate("exact location",xy=(1.05,-0.105))	
xlabel('A',size=15)
ylabel('B',size=15)
clabel(Contour,inline=True)
grid(True)

# printing the variation of error with noise in the estimate of A and B.
figure(3)
title("Q10:Variation of Error with Noise",size=15)	
xlabel('Noise standard deviation',size=15)
ylabel('Error in estimate of A and B',size=15)
plot(stdev,Err_A,label='err_A',marker='*',linestyle='dashed')
plot(stdev,Err_B,label='err_B',marker='*',linestyle='dashed')
grid(True)
legend()

# printing the variation of error with noise in log scale in the estimate of A and B.
figure(4)
title("Q11:Variation of Error with Noise on log scale",size=15)	
xlabel('stdev',size=15)
ylabel('Error in estimate of A and B',size=15)
loglog(stdev,Err_A,'bo',label='Err_A',)
errorbar(stdev, Err_A, std(Err_A), fmt='bo')
loglog(stdev,Err_B,'ro',label='Err_B')
errorbar(stdev, Err_B, std(Err_B), fmt='ro')
grid(True)
legend()

show()