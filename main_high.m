clc
clear
close all

Fs=48000;%%%采样率

t=0:1/Fs:0.02;%%时长
f=150;
w=2*pi*40;

x=sinc(w*t).*exp(1i*2*pi*f*t);     %x信号

% figure
% plot(t,abs(x));

tt=-0.01:1/Fs:0.03;

ff=650;
x_delay=sinc(w*tt).*exp(1i*2*pi*ff*tt);     %x延迟信号
% x_delay=[zeros(1,100),x,zeros(1,100)];
%x_delay=awgn(x_delay,15);
% figure
% plot(tt,abs(x_delay));


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


[a b]=max(xor_value,[],2);     %%b存放每行最大值的纵坐标，d表示了时延点数
[c d]=max(a);
d
b(d)                           %%对应的频率值 K
delta_f=Fs*b(d)/FFT_len-Fs/2   %%频差
delta_t=-d/Fs;                 %%时差
high_value= x.*conj(x).*x.*conj(x_delay(d:d+length(x)-1))-2*x_cor_conj.*x.*conj(x_delay(d:d+length(x)-1))-2*x_cor.*conj(x).*conj(x_delay(d:d+length(x)-1));
high_value= abs(high_value);
% FFT_len_2=2^14;
% value=fftshift(abs(fft(high_value,FFT_len_2)));
% [max_value, index]=max(value);
% delta_f=Fs*index/FFT_len_2-Fs/2    %直接对该行数据进行高FFT变换。



 fL=delta_f-Fs/FFT_len;
 fH=delta_f+Fs/FFT_len;
 f_search=fH-fL; f_stop=1;
 value=[];
 
 while (f_search > f_stop)
    
    i=0;
    for fk=fL:fH
        sum=0;
        for k=1:length(high_value)
            sum=sum+high_value(k)*exp(-1i*2*pi*(fk)*(k/length(high_value)+delta_t));
        end
        i=i+1;
        value(i)=abs(sum);
    end
    
    [max_value, index]=max(value);
    delta_f=fL+index-1
    fL=delta_f-Fs/FFT_len
    fH=delta_f+Fs/FFT_len
    f_search=f_search/2;
  
end

%delta_f


 
% while (f_search > f_stop)
%     temp=x.*conj(x).*x.*conj(x_delay(d))-2*x_cor_conj.*x.*conj(x_delay(d))-2*x_cor.*conj(x).*conj(x_delay(d));
%     xor_value(d,:)=fftshift(abs(fft(temp,FFT_len)));
%     
%     [a b]=max(xor_value,[],2);
%     [c d]=max(a);
%     b(d)
%     delta_f=Fs*b(d)/FFT_len-Fs/2
%     
%     f_search=f_search/2;
%     
% end

%  [aa bb]=max(xor_value,[],2);
% 
%  bb(d)
%  delta_ff=Fs*bb(d)/FFT_len-Fs/2   %%频差















