'''
Name: Hariharan P
Roll No.: EE20B042
EE2703 APL A8 DFT
'''

# Importing the necessary library
from pylab import *


N1 = 128
t1 = linspace(0,2*pi,N1+1); t1 = t1[:-1]
y1 = sin(5*t1)
Y1 = fftshift(fft(y1))/N1
w1 = linspace(-64,63,N1)

# Magnitude and phase plot for the DFT of sin(5t)
figure(0)
plot(w1,abs(Y1))
title(r"Spectrum of $\sin(5t)$")
ylabel(r"$|Y(\omega)|\rightarrow$")
xlabel(r"$\omega\rightarrow$")
xlim([-10,10])
grid(True)

figure(1)
plot(w1,angle(Y1),'ro')
ii = where(abs(Y1)>1e-3)
plot(w1[ii],angle(Y1[ii]),'go')
title(r"Phase of $\sin(5t)$")
ylabel(r"$\angle Y(\omega)\rightarrow$")
xlabel(r"$\omega\rightarrow$")
xlim([-10,10])
grid(True)


N2 = 512
t2 = linspace(-4*pi,4*pi,N2+1); t2 = t2[:-1]
y2 = (1 + 0.1*cos(t2))*cos(10*t2)
Y2 = fftshift(fft(y2))/N2
w2 = linspace(-64,64,N2+1); w2 = w2[:-1]

# Magnitude and phase plot for the DFT of (1 + 0.1cos(t))cos(10t)
figure(2)
plot(w2,abs(Y2))
title(r"Spectrum of $(1 + 0.1*cos(t))*cos(10*t)$")
ylabel(r"$|Y(\omega)|\rightarrow$")
xlabel(r"$\omega\rightarrow$")
xlim([-15,15])
grid(True)
figure(3)
plot(w2,angle(Y2),'ro')
title(r"Phase of $(1 + 0.1*cos(t))*cos(10*t)$")
ylabel(r"$\angle Y(\omega)\rightarrow$")
xlabel(r"$\omega\rightarrow$")
xlim([-15,15])
grid(True)


N3 = 512
t3 = linspace(-4*pi,4*pi,N3+1); t3 = t3[:-1]
y3 = (3*sin(t3) - sin(3*t3))/4
Y3 = fftshift(fft(y3))/N3
w3 = linspace(-64,64,N3+1); w3 = w3[:-1]

# Magnitude and phase plot for the DFT of sin^3(t)
figure(4)
plot(w3,abs(Y3))
title(r"Spectrum of $sin^{3}(t)$")
ylabel(r"$|Y(\omega)|\rightarrow$")
xlabel(r"$\omega\rightarrow$")
xlim([-15,15])
grid(True)
figure(5)
plot(w3,angle(Y3),'ro')
title(r"Phase of $sin^{3}(t)$")
ylabel(r"$\angle Y(\omega)\rightarrow$")
xlabel(r"$\omega\rightarrow$")
xlim([-15,15])
grid(True)


y4 = (3*cos(t3) + cos(3*t3))/4
Y4 = fftshift(fft(y4))/N3

# Magnitude and phase plot for the DFT of cos^3(t)
figure(6)
plot(w3,abs(Y4))
title(r"Spectrum of $cos^{3}(t)$")
ylabel(r"$|Y(\omega)|\rightarrow$")
xlabel(r"$\omega\rightarrow$")
xlim([-15,15])
grid(True)
figure(7)
plot(w3,angle(Y4),'ro')
title(r"Phase of $cos^{3}(t)$")
ylabel(r"$\angle Y(\omega)\rightarrow$")
xlabel(r"$\omega\rightarrow$")
xlim([-15,15])
grid(True)


y5 = cos(20*t3 + 5*cos(t3))
Y5 = fftshift(fft(y5))/N3

# Magnitude plot for the DFT of cos(20t + 5cos(t))
figure(8)
plot(w3,abs(Y5))
title(r"Spectrum of $cos(20t + 5cos(t))$")
ylabel(r"$|Y(\omega)|\rightarrow$")
xlabel(r"$\omega\rightarrow$")
xlim([-40,40])
grid(True)

# Phase plot for the DFT of cos(20t + 5cos(t)) ( magnitude is greater than 1e-3)
figure(9)
ii = where(abs(Y5)>=1e-3)
plot(w3[ii],angle(Y5[ii]),'go')
title(r"Phase of $cos(20t + 5cos(t))$")
ylabel(r"$\angle Y(\omega)\rightarrow$")
xlabel(r"$\omega\rightarrow$")
xlim([-40,40])
grid(True)

#The gaussian and it's precision
T = 2*pi
N = 128
Y_previous=0
tolerance=1e-6
err_t=tolerance+1

count = 0
#finding the correct window size considering tolerance using iterative loop
while err_t>tolerance:  
    x=linspace(-T/2,T/2,N+1)[:-1]
    w = linspace(-N*pi/T,N*pi/T,N+1)[:-1]
    y = exp(-0.5*x**2)
    Y=fftshift(fft(ifftshift(y)))*T/(2*pi*N)
    err_t = sum(abs(Y[::2]-Y_previous))
    Y_previous = Y
    count+=1
    T*=2
    N*=2
        
Y_expec=1/sqrt(2*pi) * exp(-w**2/2)
#computing error
error = sum(abs(abs(Y)-Y_expec))
print("True error: ",error)
print("samples = "+str(N)+"time period = "+str(T/pi) + "*pi")

figure(10)
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-5,5])
ylabel('Magnitude',size=16)
title("Estimate fft of gaussian")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
ii=where(abs(Y)>1e-3)
plot(w[ii],angle(Y[ii]),'go',lw=2)
xlim([-5,5])
ylabel("Phase",size=16)
xlabel("w",size=16)
grid(True)
show()

figure(11)
subplot(2,1,1)
plot(w,abs(Y_expec),lw=2)
xlim([-5,5])
ylabel('Magnitude',size=16)
title("True fft of gaussian")
grid(True)
subplot(2,1,2)
plot(w,angle(Y_expec),'ro',lw=2)
ii=where(abs(Y_expec)>1e-3)
plot(w[ii],angle(Y_expec[ii]),'go',lw=2)
xlim([-5,5])
ylabel("Phase",size=16)
xlabel("w",size=16)
grid(True)
show()


