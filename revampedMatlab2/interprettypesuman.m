








%load fig01_maybeok.mat;

load fig01-8888.mat

for i=1:1%change 10 to one for single config
esERR1=esERR(:,:,i);
%esgpolyERR1=esgpolyERR(:,:,:,i); %incorrect HFPOLY method TW 4/9/2020
esgpolylatestERR1 = esgpolylatestERR(:,:,:,i);
%espsERR1=espsERR(:,:,:,i); %relic of the past TW 4/9/2020
esgERR1=esgERR(:,:,i);
nsgERR1=nsgERR(i,:);
%nsgERR1=nsgERR1.';
polyERR1=polyERR(:,:,i);
psERR1=psERR(:,:,i);
esgpsERR1=esgpsERR(:,:,:,i);
fileN = ['fig',num2str(10),'_7_200.mat'];
 %save([fileN],'esERR1','esgpolyERR1','espsERR1','esgERR1','nsgERR1','polyERR1','psERR1','esgpsERR1','esgpolylatestERR1')
 save([fileN],'esERR1','esgERR1','nsgERR1','polyERR1','psERR1','esgpsERR1','esgpolylatestERR1')
%clear  esERR1 espsERR1 esgERR1 nsgERR1 polyERR1 psERR1 esgpsERR1 esgpolyERR1
clear  esERR1 esgERR1 nsgERR1 polyERR1 psERR1 esgpsERR1 esgpolyERR1
end

































