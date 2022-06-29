from pylab import *
from mpl_toolkits.mplot3d import Axes3D

def spectrum_plot(lim,n,f,t_check=0,show_check = True,t_lims = False,windowing= False,xlim1=10,title_p = r"Spectrum of $\sin\left(\sqrt{2}t\right)$",xlabel_p = r"$\omega$",ylabel1= r"$|Y|$", ylabel2 = r"Phase of $Y$"):
    if(t_lims):
        t = t_check
    else:
        t=linspace(-lim,lim,n+1)[:-1]
    dt=t[1]-t[0]
    fmax=1/dt
    y = f(t)
    if (windowing):
        m=arange(n)
        wnd=fftshift(0.54+0.46*cos(2*pi*m/n))
        y = y*wnd
    y[0]=0 
    y=fftshift(y) 
    Y=fftshift(fft(y))/float(n)
    w=linspace(-pi*fmax,pi*fmax,n+1)[:-1]
    
    mag = abs(Y)
    phase = angle(Y)
    if(show_check):
        figure()
        subplot(2,1,1)
        plot(w,mag,lw=2)
        xlim([-xlim1,xlim1])
        ylabel(ylabel1,size=16)
        title(title_p)
        grid(True)
        subplot(2,1,2)
        phase[where(mag<3e-3)] = 0
        plot(w,phase,'ro',lw=2)
        xlim([-xlim1,xlim1])
        ylabel(ylabel2,size=16)
        xlabel(xlabel_p,size=16)
        grid(True)
        show()
    return w,Y

def sinsq(t):
    return sin(sqrt(2)*t)

def cos3(t,w0=0.86):
    return (cos(w0*t))**3

def cosine(t,w0=1.5,delta=0.5):
    return cos(w0*t + delta)

def noisycosine(t,w0=1.5,delta=0.5):
    return cos(w0*t + delta) + 0.1*randn(128)

def chirp(t):
    return cos(16*(1.5 + t/(2*pi))*t)    

def omega_delta(w,Y):
    ii = where(w>0)
    omega = (sum(abs(Y[ii])**2*w[ii])/sum(abs(Y[ii])**2))
    i = abs(w-omega).argmin()
    delta = angle(Y[i])
    print ("omega = ", omega)
    print ("delta = ", delta)
   
w,Y = spectrum_plot(pi,64,sinsq,xlim1= 10,windowing=False, title_p = r"Spectrum of $\sin\left(\sqrt{2}t\right)$ without windowing ")    
w,Y = spectrum_plot(4*pi,256,sinsq,xlim1= 4,windowing=True,title_p= r"Spectrum of \sin\left(\sqrt{2}t\right)$ with windowing ")
w,Y = spectrum_plot(4*pi,256,cos3,xlim1= 4,windowing=False, title_p = r"Spectrum of $\cos^{3}(0.86t)$ without windowing ")
w,Y = spectrum_plot(4*pi,256,cos3,xlim1= 4,windowing=True, title_p = r"Spectrum of $\cos^{3}(0.86t)$ with windowing ")
w,Y = spectrum_plot(pi,128,cosine,xlim1= 4,windowing=True, title_p= r"Spectrum of $\cos(w_0t+\delta)$ with windowing ")
omega_delta(w,Y)
w,Y = spectrum_plot(pi,128,noisycosine,xlim1= 4,windowing=True, title_p = r"Spectrum of noisy $\cos(w_0t+\delta)$ with windowing ")
omega_delta(w,Y)
w,Y= spectrum_plot(pi,1024,chirp,xlim1= 60,windowing=False, title_p = r"Spectrum of chirp function without windowing")
w,Y = spectrum_plot(pi,1024,chirp,xlim1= 60,windowing=True,title_p = r"Spectrum of chirp function with windowing")

t = linspace(-pi,pi,1025); t = t[:-1]
dt = t[1]-t[0]; fmax = 1/dt
t_array = split(t,16)
Y_mag = zeros((16,64))
Y_phase = zeros((16,64))
for i in range(len(t_array)):
    w,Y = spectrum_plot(0,64,f = chirp,show_check = False,t_lims=True,t_check = t_array[i],xlim1= 60,windowing=True, title_p = r"Spectrum of chirp function")
    Y_mag[i] =  abs(Y)
    Y_phase[i] =  angle(Y)

t= t[::64]	
w = linspace(-fmax*pi,fmax*pi,65); w = w[:-1]
t,w = meshgrid(t,w)

fig1 = figure(9)
ax = fig1.add_subplot(111, projection='3d')
surf=ax.plot_surface(w,t,Y_mag.T,cmap='viridis',linewidth=0, antialiased=False)
fig1.colorbar(surf, shrink=0.5, aspect=5)
ax.set_title('magnitude surface plot');
ylabel(r"$\omega\rightarrow$")
xlabel(r"$t\rightarrow$")
show()
fig1 = figure(10)
ax = fig1.add_subplot(111, projection='3d')
surf=ax.plot_surface(w,t,Y_phase.T,cmap='viridis',linewidth=0, antialiased=False)
fig1.colorbar(surf, shrink=0.5, aspect=5)
ax.set_title('phase surface plot');
ylabel(r"$\omega\rightarrow$")
xlabel(r"$t\rightarrow$")
show()