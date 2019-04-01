clc
clear
close all

Fs=48000;%%%采样率

t=0:1/Fs:0.02;%%时长
f=150;
w=2*pi*40;

x=sinc(w*t).*exp(1i*2*pi*f*t);     %x信号
x=awgn(x,7);
 figure
 plot(t,abs(x));

tt=-0.01:1/Fs:0.03;

ff=650;
x_delay=sinc(w*tt).*exp(1i*2*pi*ff*tt);     %x延迟信号
% x_delay=[zeros(1,100),x,zeros(1,100)];
x_delay=awgn(x_delay,7);
 figure
 plot(tt,abs(x_delay));


%%%时差为0.01 频差为-500
xor_len=length(x_delay)-length(x)+1;
FFT_len=2^12;

xor_value=zeros(xor_len,FFT_len);

for  ii=1:length(x)
    x_cor=0;
    x_cor=x_cor+x(ii)*x(ii); 
end
x_cor=x_cor/length(x);

for  jj=1:length(conj(x))
    x_cor_conj=0;
    x_cor_conj=x_cor_conj+conj(x(jj))*conj(x(jj)); 
end
x_cor_conj=x_cor_conj/length(conj(x)); 

for i=1:xor_len
    temp=x.*conj(x).*x.*conj(x_delay(i:i+length(x)-1))-2*x_cor_conj.*x.*conj(x_delay(i:i+length(x)-1))-2*x_cor.*conj(x).*conj(x_delay(i:i+length(x)-1));
    
    %temp=x.*conj(x_delay(i:i+length(x)-1));
    xor_value(i,:)=fftshift(abs(fft(temp,FFT_len))); 
end

[a b]=max(xor_value,[],2);
[c d]=max(a);                   %%d为时延的点数，相当于上面的i

b(d)
delta_f=Fs*b(d)/FFT_len-Fs/2    %%频差
delta_t=-d/Fs                   %%时差

figure
% TT=-0.01:1/Fs:0.01;
TT=-0.02:1/Fs:0.00;
FF=(0:Fs/FFT_len:Fs-1)-Fs/2;
[X,Y]=meshgrid(FF,TT);
mesh(X,Y,xor_value)
xlabel('f');
ylabel('t');
title('CAF')

 


















