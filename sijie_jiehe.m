
clear
close all

Fs=256;   %%������
N=128;
n=N/2:3/2*N-1;%fc=100;w=2*pi*fc;  %�ز���fc
t=n/Fs;     %%ʱ��
f1=61;

x=exp(1i*(2*pi*f1*t+(2*rand(1,1)-1)*2*pi));     %x�ź�
 x=awgn(x,5);
h=Blackman(N); 
h1=Blackman(2*N);
% for m=1:N
%     x(m)=x(m)*h(m);
% end
nn=0:2*N-1;
tt=nn/Fs;

f2=20;
x_delay=exp(1i*(2*pi*f2*tt+(2*rand(1,1)-1)*2*pi));     %x�ӳ��ź�
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
[c d]=max(a);                  %%dΪʱ�ӵĵ������൱�������i

b(d)
delta_f=(b(d)-1)*Fs/N-Fs/2    %%Ƶ��/
delta_t=-d/Fs   

y=x.*conj(x_delay(d:d+length(x)-1));
y1=y(1:N/2);
y2=y(N/2+1:N);
y1fft=fft(y1,N/2);
y2fft=fft(y2,N/2);
[a1,b1]=max(y1fft);
[a2,b2]=max(y2fft);

f1max=(max(b1,b2)-1)*Fs/(N/2)     %������ߴֲ�Ƶ��  (b1-1)
ph1=angle(y1fft(1:N/2));    
ph2=angle(y2fft(1:N/2));
b1
b2

angle11=ph1(min(b1,b2))            %y1�ź���������ߴ�����λ
angle22=ph2(min(b1,b2))
angle33=angle22-angle11      %%�ó������ź���DFT������ߴ�����λ��
f_err=angle33/(2*pi)*2*Fs/N   %��ff��������߹��Ƶ�Ƶ�ʵ�ƫ����й���
f_est=f_err+f1max          %���Ƶ�Ƶ��


T=(N)/Fs;           %  T�Ĺ�������N����
fe=f_err/(2*Fs/N)    %���Ƶƫ����������0.5֮��  angle33/2pi
snr=(0:1:30);
% snr=10.^(snr./10);
%Ƶ�ʹ��ƾ������������ֵ
ro_error=1./(((N/2)*10.^(snr./10)).^(1/2))*1/(abs(sinc(fe))*pi*T);  

%plot(snr,ro_error,'b-'); 
%semilogy(snr,ro_error,'k--');

axis([0 30 -2 10]);
% legend('Theoretical result','Simulation result','CR Lower limit');
%axis([-inf inf ]) %xmin��x��С��xmax��x���ymin��ymax����








