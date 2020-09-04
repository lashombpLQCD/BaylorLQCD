clear all;
echo off

fileF = 'hopping-matrix.LOG';
%fileM = (['/data/barals/qqcd/scratchhoppingmatrix/hopteninone.mat']);

  fid=fopen([fileF],'r');

% fid=fopen(fileF,'r');
el = 12939264;
totalconfig=1;
  numelem = el*totalconfig;
  matrixsize = 12 * 12^3 *16;
  config =1;
  % data holders
  ii   = zeros(numelem,totalconfig);
  jj   = zeros(numelem,totalconfig);
  elem = zeros(numelem,totalconfig);
  currentline =0;
 % for j=(numelem*(x-1)+1):(numelem*x)
 % for j=(numelem*(x-1)+1):(numelem*x)
  for j=1:numelem
currentline=currentline+1;
%    if (mod(j,20000) == 0)
 %       disp(['j: ',num2str(j)])
  %  end

      dat=fscanf(fid,'%i %i %f %f',[1 4]);
% p = (j-numelem*(x-1));
      ii(j) = dat(1);
      jj(j) = dat(2);
      elem(j) = dat(3) + 1i*dat(4);
if (ii(j)<=0)
      disp(['Currentline in ii not positive: ',int2str(currentline)])
elseif (jj(j)<=0)
      disp(['Currentline in jj not positive: ',int2str(currentline)])
end 



 if (mod(j,el) == 0)
      disp(['New configuration starts here: ',num2str(j)])
      config = config +1;




    end


  end%for j 


oneconfig = numelem/(config-1)
config=1;

for (k=oneconfig:oneconfig:numelem)
    for j =1:oneconfig
      aa(j,config)=ii(j+k-oneconfig);  
     bb(j,config)=jj(j+k-oneconfig);  
      cc(j,config)=elem(j+k-oneconfig);  
    end
        config =config+1;
end


%  C=sparse(matrixsize,matrixsize,config);
for config =1:totalconfig
if (config<10)
fileM = (['h0',num2str(config),'.mat']);
else
fileM = (['h',num2str(config),'.mat']);
end
  B = sparse(aa(:,config),bb(:,config),cc(:,config),matrixsize,matrixsize);





%Change color and dirac in hopping matrices. To speed up transpose is used. nothing fancy.


%old/new: 1/1 2/5 3/9 4/2 5/6 6/10 7/3 8/7 9/11 10/4 11/8 12/12
% looks like this changes color as leading index to dirac as leading index - TW 4/16/2020


for j=1:12:matrixsize

% if (mod(j,2000) == 0)
      disp(['j2: ',num2str(j)])
 %   end

tempcol1=B(:,j+1);
tempcol2=B(:,j+2);
tempcol3=B(:,j+4);
tempcol4=B(:,j+5);
tempcol5=B(:,j+8);

B(:,j+1)=B(:,j+3);
B(:,j+2)=B(:,j+6);
B(:,j+3)=B(:,j+9);
B(:,j+4)=tempcol1;
B(:,j+5)=tempcol3;
B(:,j+6)=B(:,j+7);
B(:,j+7)=B(:,j+10);
B(:,j+8)=tempcol2;
B(:,j+9)=tempcol4;
B(:,j+10)=tempcol5;
end
disp('after for loop')
%old/new: 1/1 2/5 3/9 4/2 5/6 6/10 7/3 8/7 9/11 10/4 11/8 12/12

  disp(['j3: just started first transpose'])
tic;
B=B.';
toc;
  disp(['j3: just did first transpose'])
tic;
for j=1:12:matrixsize
  disp(['j3: ',num2str(j)])

tempcol1=B(:,j+1);
tempcol2=B(:,j+2);
tempcol3=B(:,j+4);
tempcol4=B(:,j+5);
tempcol5=B(:,j+8);

B(:,j+1)=B(:,j+3);
B(:,j+2)=B(:,j+6);
B(:,j+3)=B(:,j+9);
B(:,j+4)=tempcol1;
B(:,j+5)=tempcol3;
B(:,j+6)=B(:,j+7);
B(:,j+7)=B(:,j+10);
B(:,j+8)=tempcol2;
B(:,j+9)=tempcol4;
B(:,j+10)=tempcol5;
end
toc;
 disp(['for final  transpose'])
tic;
B =B.';
toc;



  save([fileM],'B');
 
clear B;
end


 % fclose(fid);
  clear B p pp ans dat fid j numelem matrixsize ii jj elem tempcol1 tempcol2 tempcol3 tempcol4 tempcol5 temprow1 temprow2 temprow3 temprow4 temprow5

