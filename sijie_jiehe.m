
clear
close all

Fs=256;   %%采样率
N=128;
n=N/2:3/2*N-1;%fc=100;w=2*pi*fc;  %载波是fc
t=n/Fs;     %%时长
f1=61;

x=exp(1i*(2*pi*f1*t+(2*rand(1,1)-1)*2*pi));     %x信号
 x=awgn(x,5);
h=Blackman(N); 
h1=Blackman(2*N);
% for m=1:N
%     x(m)=x(m)*h(m);
% end
nn=0:2*N-1;
tt=nn/Fs;

f2=20;
x_delay=exp(1i*(2*pi*f2*tt+(2*rand(1,1)-1)*2*pi));     %x延迟信号
 x_delay=awgn(x_delay,0);
% for m=1:2*N
%     x_delay(m)=x_delay(m)*h1(m);
% end


xor_len=length(x_delay)-length(x);
xor_value=zeros(xor_len,N);

x_cor=0;
for  ii=1:length(x)
   
    x_cor=x_cor+x(ii)*x(ii); 
end
x_cor=x_cor/length(x);

x_cor_conj=0;
for  jj=1:length(conj(x))
    
    x_cor_conj=x_cor_conj+conj(x(jj))*conj(x(jj)); 
end
x_cor_conj=x_cor_conj/length(conj(x)); 

for i=1:xor_len
    temp=x.*conj(x).*x.*conj(x_delay(i:i+length(x)-1))-2*x_cor_conj.*x.*conj(x_delay(i:i+length(x)-1))-x_cor.*conj(x).*conj(x_delay(i:i+length(x)-1));
    
    %temp=x.*conj(x_delay(i:i+length(x)-1));
    xor_value(i,:)=fftshift(abs(fft(temp,N))); 
end

[a b]=max(xor_value,[],2);
[c d]=max(a);                  %%d为时延的点数，相当于上面的i

b(d)
delta_f=(b(d)-1)*Fs/N-Fs/2    %%频差/
delta_t=-d/Fs   

y=x.*conj(x_delay(d:d+length(x)-1));
y1=y(1:N/2);
y2=y(N/2+1:N);
y1fft=fft(y1,N/2);
y2fft=fft(y2,N/2);
[a1,b1]=max(y1fft);
[a2,b2]=max(y2fft);

f1max=(max(b1,b2)-1)*Fs/(N/2)     %最大谱线粗测频率  (b1-1)
ph1=angle(y1fft(1:N/2));    
ph2=angle(y2fft(1:N/2));
b1
b2

angle11=ph1(min(b1,b2))            %y1信号在最大谱线处的相位
angle22=ph2(min(b1,b2))
angle33=angle22-angle11      %%得出两段信号在DFT最大谱线处的相位差
f_err=angle33/(2*pi)*2*Fs/N   %对ff与最大谱线估计的频率的偏差进行估计
f_est=f_err+f1max          %估计的频率


T=(N)/Fs;           %  T的估计整数N个点
fe=f_err/(2*Fs/N)    %相对频偏，介于正负0.5之间  angle33/2pi
snr=(0:1:30);
% snr=10.^(snr./10);
%频率估计均方根误差理论值
ro_error=1./(((N/2)*10.^(snr./10)).^(1/2))*1/(abs(sinc(fe))*pi*T);  

%plot(snr,ro_error,'b-'); 
%semilogy(snr,ro_error,'k--');

axis([0 30 -2 10]);
% legend('Theoretical result','Simulation result','CR Lower limit');
%axis([-inf inf ]) %xmin是x最小，xmax是x最大，ymin，ymax类似








