
function f=sijie_hanshu()
clc
clear
close all

Fs=48000;%%%������
t=0:1/Fs:0.02;%%ʱ��
f=1200;
w=2*pi*40;

x=5*sinc(w*t).*exp(1i*2*pi*f*t);     %x�ź�
x=awgn(x,-5);
% figure
% plot(t,abs(x));

tt=-0.01:1/Fs:0.03;
ff=200;
x_delay=sinc(w*tt).*exp(1i*2*pi*ff*tt);     %x�ź�
x_delay=awgn(x_delay,-5);
% figure
% plot(tt,abs(x_delay));


%%%ʱ��Ϊ0.01 Ƶ��Ϊ1000
xor_len=length(x_delay)-length(x)+1;
FFT_len=2^11;

xor_value=zeros(xor_len,FFT_len);     %%���й�ϵ
%%%��������ά��ģ��
for i=1:xor_len
    temp=x.*conj(x_delay(i:i+length(x)-1));
    xor_value(i,:)=fftshift(abs(fft(temp,FFT_len)));  %%���Ժ���exp��1i*2*pi*(Fd-f)*n����
end

[a b]=max(xor_value,[],2);  %%ÿ�����ֵ���a�� b
[c d]=max(a);

b(d);                      %%
delta_f=Fs*b(d)/FFT_len-Fs/2;  %%Ƶ��
delta_t=d/Fs;                 %%ʱ��

f=delta_f


