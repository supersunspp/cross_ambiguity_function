clc
clear
close all

Fs=48000;%%%������

t=0:1/Fs:0.02;%%ʱ��
f=1500;
w=2*pi*40;

x=sinc(w*t).*exp(1i*2*pi*f*t);     %x�ź�
x=awgn(x,20);
figure
plot(t,abs(x));

tt=-0.01:1/Fs:0.03;

ff=2000;
x_delay=sinc(w*tt).*exp(1i*2*pi*ff*tt);     %x�ź�

x_delay=awgn(x_delay,20);
figure
plot(tt,abs(x_delay));


%%%ʱ��Ϊ0.01 Ƶ��Ϊ-500
xor_len=length(x_delay)-length(x)+1;

FFT_len=2^12;

xor_value=zeros(xor_len,FFT_len);
%%%��������ά��ģ��
for i=1:xor_len
    temp=x.*conj(x_delay(i:i+length(x)-1));
    xor_value(i,:)=fftshift(abs(fft(temp,FFT_len))); 
end

[a b]=max(xor_value,[],2);
[c d]=max(a);

b(d)
delta_f=Fs*b(d)/FFT_len-Fs/2    %%Ƶ��
delta_t=-d/Fs                   %%ʱ��

figure
TT=-0.01:1/Fs:0.01;
TT=-0.02:1/Fs:0.00;
FF=(0:Fs/FFT_len:Fs-1)-Fs/2;
[X,Y]=meshgrid(FF,TT);
mesh(X,Y,xor_value)
xlabel('f');
ylabel('t');
title('CAF')
