config=10
rnk = 331776;
numlines=rnk*config
tic;
  for i=1:200
      M = dlmread(['noise-',num2str(i),'.LOG']);
      B(:,:,i)=M(:,2:3);
    clear M   
i
  end

  for k = 1:200 
     for j=1:numlines
      C(j,k)=B(j,1,k)+1i*B(j,2,k);
     end  
  end
    clear B

  for i=1:config
    fileM = (['z2--',num2str(i),'.mat']);
    num = rnk*(i-1)
    z2(1:rnk,:) = C((num+1:rnk+num),:);

  z2=z2.';


for j=1:12:rnk
  disp(['j3: ',num2str(j)])

  tempcol1=z2(:,j+1);
  tempcol2=z2(:,j+2);
  tempcol3=z2(:,j+4);
  tempcol4=z2(:,j+5);
  tempcol5=z2(:,j+8);

  z2(:,j+1)=z2(:,j+3);
  z2(:,j+2)=z2(:,j+6);
  z2(:,j+3)=z2(:,j+9);
  z2(:,j+4)=tempcol1;
  z2(:,j+5)=tempcol3;
  z2(:,j+6)=z2(:,j+7);
  z2(:,j+7)=z2(:,j+10);
  z2(:,j+8)=tempcol2;
  z2(:,j+9)=tempcol4;
  z2(:,j+10)=tempcol5;
 end%j

   z2 =z2.';


    save([fileM],'z2');
   clear z2;
  end%config


clear C B




