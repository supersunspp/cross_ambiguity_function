cha=zeros(1,100);
for  v=1:100
    f=sijie_hanshu();
   cha(v)=abs(1000-f);
end
aver_value=mean(cha)   %��ֵ
ro_error2=(sum(cha.^2)/100)^(1/2)  %���������ֵ