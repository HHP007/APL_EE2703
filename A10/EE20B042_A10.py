import numpy as np
import scipy.signal as sig 
from pylab import *
import csv
def plot_fun(x,y,x1='w',y1='Magnitude',t1='xyz'):
    plot(x,y)
    xlabel(x1)
    ylabel(y1)
    title(t1)
    grid(True)
    show()

details =np.zeros(12)
i = 0
with open("h.csv", 'r') as file:  
    for line in file:
        details[i] = float(line)
        i+=1

w,h = sig.freqz(details)
plot_fun(w,abs(h),'w','Magnitude','Low pass filter')
plot_fun(w,angle(h),'w','Phase','Low pass filter')

n = arange(2**10)
x = cos(0.2*pi*n) + cos(0.85*pi*n)
plot_fun(n,x,'n','Amplitude','Input signal')

y = np.zeros(len(x))
for i in arange(len(x)):
    for k in arange(len(details)):
        y[i]+=x[i-k]*details[k]        
plot_fun(n,y,'n','Amplitude','Output of linear convolution')

y=ifft(fft(x)*fft(concatenate((details,zeros(len(x)-len(details))))))
plot_fun(n,real(y),'n','Amplitude','Output of circular convolution')

def circular_conv(x,h):
    P = len(h)
    n_ = int(ceil(log2(P)))
    h_ = np.concatenate((h,np.zeros(int(2**n_)-P)))
    P = len(h_)
    n1 = int(ceil(len(x)/2**n_))
    x_ = np.concatenate((x,np.zeros(n1*(int(2**n_))-len(x))))
    y = np.zeros(len(x_)+len(h_)-1)
    for i in range(n1):
        temp = np.concatenate((x_[i*P:(i+1)*P],np.zeros(P-1)))
        y[i*P:(i+1)*P+P-1] += np.fft.ifft(np.fft.fft(temp) * np.fft.fft( np.concatenate( (h_,np.zeros(len(temp)-len(h_))) ))).real
    return y

y = circular_conv(x,details)
plot_fun(n,real(y[:1024]),'n','Amplitude','Output of linear convolution using circular convolution')

lines = []
with open("x1.csv",'r') as file2:
    csvreader = csv.reader(file2)
    for row in csvreader:
        lines.append(row)
lines2 = []
for line in lines:
    line = list(line[0])
    try :
        line[line.index('i')]='j'
        lines2.append(line)
    except ValueError:
        lines2.append(line)
        continue
x = [complex(''.join(line)) for line in lines2]
x2 = np.roll(x,5)
cor = np.fft.ifftshift(np.correlate(x2,x,'full'))
xlim(0,25)
plot_fun(linspace(0,len(cor)-1,len(cor)),abs(cor),'t','Correlation','auto-correlation')