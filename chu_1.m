clc
clear all
%��λ�ڲ巨
Fs=1024;   %%������
N=2^9;
n=0:N-1;
t=n/Fs;     %%ʱ��
angle1=pi/2; angle2=pi/4;
f1=341; 
%noise=sqrt(1)*randn(1,N)+1i*sqrt(1)*randn(1,N);

x=exp(1i*2*pi*f1*t+angle1);
x=awgn(x,0);

f2=300;
nn=N/2:3/2*N-1;
tt=nn/Fs;
x_delay=exp(1i*2*pi*f2*tt+angle2);
x_delay=awgn(x_delay,0);

y=x./x_delay;
y1=y(1:N/2);
y2=y(N/2+1:N);
y1fft=fft(y1,N/2);
y2fft=fft(y2,N/2);
[a1,b1]=max(y1fft);
f1max=(b1-1)*Fs/(N/2)   %������ߴֲ���Ƶ��
[a2,b2]=max(y2fft);
ph1=angle(y1fft(1:N/2));    
ph2=angle(y2fft(1:N/2));
b1
b2

angle11=ph1(b1);            %y1�ź���������ߴ�����λ
angle22=ph2(b2);
angle=angle22-angle11;      %%�ó������ź���DFT������ߴ�����λ��
f_err=angle/(2*pi)*2*Fs/N   %��ff��������߹��Ƶ�Ƶ�ʵ�ƫ����й���
f_est=f_err+f1max           %���Ƶ�Ƶ��














