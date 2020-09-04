iopPLOT = 5;
korder  = 7;
config  =10;
nop = 5; %5 for matlab sub prog, 9 for fortran sub prog
%korder = korder + 1;
neig =0;
neigH =0;
% loag raw data
 noises = 200;

% fig01 = load(['C:\Users\Suman\Desktop\interpretfiles\1560fortranlarge\figsmooth1.mat']);
% fig02 = load(['C:\Users\Suman\Desktop\interpretfiles\1560fortranlarge\figsmooth2.mat']);
% fig03 = load(['C:\Users\Suman\Desktop\interpretfiles\1560fortranlarge\figsmooth3.mat']);
% fig04 = load(['C:\Users\Suman\Desktop\interpretfiles\1560fortranlarge\figsmooth4.mat']);
% fig05 = load(['C:\Users\Suman\Desktop\interpretfiles\1560fortranlarge\figsmooth5.mat']);
% fig06 = load(['C:\Users\Suman\Desktop\interpretfiles\1560fortranlarge\figsmooth6.mat']);
% fig07 = load(['C:\Users\Suman\Desktop\interpretfiles\1560fortranlarge\figsmooth7.mat']);
% fig08 = load(['C:\Users\Suman\Desktop\interpretfiles\1560fortranlarge\figsmooth8.mat']);
% fig09 = load(['C:\Users\Suman\Desktop\interpretfiles\1560fortranlarge\figsmooth9.mat']);
% fig10 = load(['C:\Users\Suman\Desktop\interpretfiles\1560fortranlarge\figsmooth10.mat']);
%fig01 = load(['C:\Users\Admin\Desktop\For Travis\fig1_maybeok.mat']);
%fig01 = load(['C:\Users\Travis_Whyte\Desktop\For Travis\fig1_8_100.mat']);
fig01 = load(['C:\Users\Travis_Whyte\Desktop\For Travis\fig1_7_200.mat']);
fig02 = load(['C:\Users\Travis_Whyte\Desktop\For Travis\fig2_7_200.mat']);
fig03 = load(['C:\Users\Travis_Whyte\Desktop\For Travis\fig3_7_200.mat']);
fig04 = load(['C:\Users\Travis_Whyte\Desktop\For Travis\fig4_7_200.mat']);
fig05 = load(['C:\Users\Travis_Whyte\Desktop\For Travis\fig5_7_200.mat']);
fig06 = load(['C:\Users\Travis_Whyte\Desktop\For Travis\fig6_7_200.mat']);
fig07 = load(['C:\Users\Travis_Whyte\Desktop\For Travis\fig7_7_200.mat']);
fig08 = load(['C:\Users\Travis_Whyte\Desktop\For Travis\fig8_7_200.mat']);
fig09 = load(['C:\Users\Travis_Whyte\Desktop\For Travis\fig9_7_200.mat']);
fig10 = load(['C:\Users\Travis_Whyte\Desktop\For Travis\fig10_7_200.mat']);
%neig = zeros(143,1);
%neig  = fig022.neig;
neig =132;
%neigH = zeros(143,1);
%neigH(01) = fig01.neigH;
%neigH = fig022.neigH;
neigH =132;
%neigH(02) = fig02.neigH;
%neigH(03) = fig03.neigH;
%neigH(04) = fig04.neigH;
%neigH(05) = fig05.neigH;
%neigH(06) = fig06.neigH;
%neigH(07) = fig07.neigH;
%neigH(08) = fig08.neigH;
%neigH(09) = fig09.neigH;
%neigH(10) = fig10.neigH;
%neigH = min(neigH);

nmax = max(neig,neigH);

% average ESg aka HFES
esgERRAVG = zeros(neigH,nop);
for iop=1:nop
  for ieig=1:neigH
    val = 0;
    val = val + fig04.esgERR1(ieig,iop);
    % Actually esgERRAVG contains the same element as esgERR
     val = val + fig02.esgERR1(ieig,iop);
     val = val + fig03.esgERR1(ieig,iop);
    val = val + fig01.esgERR1(ieig,iop);
     val = val + fig05.esgERR1(ieig,iop);
     val = val + fig06.esgERR1(ieig,iop);
     val = val + fig07.esgERR1(ieig,iop);
     val = val + fig08.esgERR1(ieig,iop);
     val = val + fig09.esgERR1(ieig,iop);
     val = val + fig10.esgERR1(ieig,iop);
    %esgERRAVG(ieig,iop) = val/10;
    esgERRAVG(ieig,iop) = val/config;
    % For number of trials  =1 esgERR is same as esgERRAVG
  end
end

% average ES
esERRAVG = zeros(neig,nop);
for iop=1:nop
  for ieig=1:neig
    val = 0;
    val = val + fig04.esERR1(ieig,iop);
     val = val + fig02.esERR1(ieig,iop);
     val = val + fig03.esERR1(ieig,iop);
    val = val + fig01.esERR1(ieig,iop);
     val = val + fig05.esERR1(ieig,iop);
     val = val + fig06.esERR1(ieig,iop);
     val = val + fig07.esERR1(ieig,iop);
     val = val + fig08.esERR1(ieig,iop);
     val = val + fig09.esERR1(ieig,iop);
     val = val + fig10.esERR1(ieig,iop);
    %esERRAVG(ieig,iop) = val/10;
    esERRAVG(ieig,iop) = val/config;
  end
end

% average NS
nsgERRAVG = zeros(nmax,nop);
for iop=1:nop
  val = 0;
  val = val + fig04.nsgERR1(1,iop);
  val = val + fig01.nsgERR1(1,iop);
   val = val + fig03.nsgERR1(1,iop);
   val = val + fig02.nsgERR1(1,iop);
   val = val + fig05.nsgERR1(1,iop);
   val = val + fig06.nsgERR1(1,iop);
   val = val + fig07.nsgERR1(1,iop);
   val = val + fig08.nsgERR1(1,iop);
   val = val + fig09.nsgERR1(1,iop);
   val = val + fig10.nsgERR1(1,iop);
  %nsgERRAVG(:,iop) = val/10 * ones(nmax,1);
  nsgERRAVG(:,iop) = val/config * ones(nmax,1);
  % Basically each column of wnsgERRAVG is same as nsgERR 
end

% average PS
psERRAVG = zeros(nmax,korder,nop);
for ksub = 1:korder
  for iop=1:nop
    val = 0;
    %val = val + fig022.psERR1(ksub,iop);
     val = val + fig02.psERR1(ksub,iop);
     val = val + fig03.psERR1(ksub,iop);
    val = val + fig01.psERR1(ksub,iop);
     val = val + fig04.psERR1(ksub,iop);
     val = val + fig05.psERR1(ksub,iop);
     val = val + fig06.psERR1(ksub,iop);
     val = val + fig07.psERR1(ksub,iop);
     val = val + fig08.psERR1(ksub,iop);
     val = val + fig09.psERR1(ksub,iop);
     val = val + fig10.psERR1(ksub,iop);
   %psERRAVG(:,ksub,iop) = val/10 * ones(nmax,1);
    psERRAVG(:,ksub,iop) = val/config * ones(nmax,1);
    % 5 such matrix named psERRAVG where a column has same elments 
  end
end
  A=psERRAVG(neig,korder,iopPLOT)

% average poly
 polyERRAVG = zeros (nmax,korder,nop);
 for ksub = 1:korder
     for iop =1:nop
         val = 0;
         val = val + fig01.polyERR1(ksub,iop);
          val = val + fig02.polyERR1(ksub,iop);
          val = val + fig03.polyERR1(ksub,iop);
          val = val + fig04.polyERR1(ksub,iop);
          val = val + fig05.polyERR1(ksub,iop);
          val = val + fig06.polyERR1(ksub,iop);
          val = val + fig07.polyERR1(ksub,iop);
          val = val + fig08.polyERR1(ksub,iop);
          val = val + fig09.polyERR1(ksub,iop);
          val = val + fig10.polyERR1(ksub,iop);
         polyERRAVG(:,ksub,iop)= val/config * ones(nmax,1);
     end 
 end


% average ESPS 
% relic of the past, no longer needed TW 4/9/2020
%*********************************************
%espsERRAVG = zeros(neig,korder,nop);
%for ksub = 1:korder
%  for iop=1:nop
%    for ieig=1:neig
%      val = 0;
%       val = val + fig02.espsERR1(ieig,ksub,iop);
%      val = val + fig01.espsERR1(ieig,ksub,iop);
%       val = val + fig03.espsERR1(ieig,ksub,iop);
%       val = val + fig04.espsERR1(ieig,ksub,iop);
%       val = val + fig05.espsERR1(ieig,ksub,iop);
%       val = val + fig06.espsERR1(ieig,ksub,iop);
%       val = val + fig07.espsERR1(ieig,ksub,iop);
%       val = val + fig08.espsERR1(ieig,ksub,iop);
%       val = val + fig09.espsERR1(ieig,ksub,iop);
%       val = val + fig10.espsERR1(ieig,ksub,iop);
      
     % espsERRAVG(ieig,ksub,iop) = val/10;
%      espsERRAVG(ieig,ksub,iop) = val/config;
      % espsERRAVG is same as espsERR for number of trials equal to one
%    end
%  end
%end
%*************************************************

% average ESgPS aka HFPS
esgpsERRAVG = zeros(neigH,korder,nop);
for ksub = 1:korder
  for iop=1:nop
    for ieig=1:neigH
      val = 0;
      val = val + fig01.esgpsERR1(ieig,ksub,iop);
       val = val + fig02.esgpsERR1(ieig,ksub,iop);
       val = val + fig03.esgpsERR1(ieig,ksub,iop);
       val = val + fig04.esgpsERR1(ieig,ksub,iop);
       val = val + fig05.esgpsERR1(ieig,ksub,iop);
       val = val + fig06.esgpsERR1(ieig,ksub,iop);
       val = val + fig07.esgpsERR1(ieig,ksub,iop);
       val = val + fig08.esgpsERR1(ieig,ksub,iop);
       val = val + fig09.esgpsERR1(ieig,ksub,iop);
       val = val + fig10.esgpsERR1(ieig,ksub,iop);
      %esgpsERRAVG(ieig,ksub,iop) = val/10;
      esgpsERRAVG(ieig,ksub,iop) = val/config;
      % esgpsERRAVG is same as esgpsERR for number of trials equal 1
    end
  end
end

%average esgpoly
%no longer needed, this was incorrect TW 4/9/2020
%*************************************************
%esgpolyERRAVG = zeros(neigH,korder,nop);
%for ksub = 1:korder
%  for iop=1:nop
%    for ieig=1:neigH
%      val = 0;
%       val = val + fig03.esgpolyERR1(ieig,ksub,iop);
%      val = val + fig01.esgpolyERR1(ieig,ksub,iop);
%       val = val + fig02.esgpolyERR1(ieig,ksub,iop);
%       val = val + fig04.esgpolyERR1(ieig,ksub,iop);
%       val = val + fig05.esgpolyERR1(ieig,ksub,iop);
%       val = val + fig06.esgpolyERR1(ieig,ksub,iop);
%       val = val + fig07.esgpolyERR1(ieig,ksub,iop);
%       val = val + fig08.esgpolyERR1(ieig,ksub,iop);
%       val = val + fig09.esgpolyERR1(ieig,ksub,iop);
%       val = val + fig10.esgpolyERR1(ieig,ksub,iop);
%      esgpolyERRAVG(ieig,ksub,iop) = val/config;
     
%    end
%  end
%end
% B=esgpolyERRAVG(neig,korder,iopPLOT)
%************************************************

% -- HFPOLY --
esgpolylatestERRAVG = zeros(neigH,korder,nop);
for ksub = 1:korder
  for iop=1:nop
    for ieig=1:neigH
      val = 0;
       val = val + fig03.esgpolyERR1(ieig,ksub,iop);
      val = val + fig01.esgpolylatestERR1(ieig,ksub,iop);
       val = val + fig02.esgpolyERR1(ieig,ksub,iop);
       val = val + fig04.esgpolyERR1(ieig,ksub,iop);
       val = val + fig05.esgpolyERR1(ieig,ksub,iop);
       val = val + fig06.esgpolyERR1(ieig,ksub,iop);
       val = val + fig07.esgpolyERR1(ieig,ksub,iop);
       val = val + fig08.esgpolyERR1(ieig,ksub,iop);
       val = val + fig09.esgpolyERR1(ieig,ksub,iop);
       val = val + fig10.esgpolyERR1(ieig,ksub,iop);
      esgpolylatestERRAVG(ieig,ksub,iop) = val/config;
     
    end
  end
end
B=esgpolylatestERRAVG(neig,korder,iopPLOT)

xmax = 1:nmax;
xeig = 1:neig;
xeigH = 1:neigH;

if (1==1)
  fig = figure;
  hold on


    plot(xmax,nsgERRAVG(:,iopPLOT),'-hb','LineWidth',1);
% It should be a straight line
  plot(xeig,esERRAVG(:,iopPLOT),'-og', 'LineWidth',1);
  %It has different elements derived from esERR
  plot(xeigH,esgERRAVG(:,iopPLOT),'-+m','LineWidth',1);

  plot(xmax,psERRAVG(:,korder,iopPLOT),'-*r','LineWidth',1);
  plot(xmax,polyERRAVG(:,korder,iopPLOT),'-*g','LineWidth',1);

 %This line just picks one of possible twenty five columns(by changing korder and iopPLOT)
  plot(xeigH,esgpsERRAVG(:,korder,iopPLOT),'-dk','LineWidth',1);
  %plot(xeigH,esgpolyERRAVG(:,korder,iopPLOT),'-dm','LineWidth',1);
  plot(xeigH,esgpolylatestERRAVG(:,korder,iopPLOT),'-dm','LineWidth',1);
  % plot(xeigH,espsERRAVG(:,korder,iopPLOT),'-dr','LineWidth',1); 
 
set(gca,'yscale','log')
  % legend('NS','ES','HFES','PS','POLY','HFPS combo','HFPOLY combo');
%   legend('NS','ES','HFES','PS','POLY','HFPS combo','ESPS combo');
  legend('NS','ES','HFES','PS','POLY','HFPS combo','HFPOLY');
  legend('Location','Best');
%  legend('NS','ES','HFES','PS','HFPS combo','ESPS combo');
%title(['J_{' num2str(iopPLOT) '}, ' num2str(korder) 'th Order, ' num2str(noises) ' noises'])
if (iopPLOT == 5)
    title(['Scalar Operator, ' num2str(korder) 'th Order, ' num2str(noises) ' noises'])
    ylim([4.5*10^(-4),1.8*10^(-3)]);
else
    title(['J_{' num2str(iopPLOT) '} Operator, ' num2str(korder) 'th Order, ' num2str(noises) ' noises'])
    ylim([2.5*10^(-4),1.8*10^(-3)]);
end

 % if (iopPLOT<5)
   % title([num2str(latsize),'  bc',bc,'  kappa: 0.',kappa,'   operator: \gamma_',num2str(iopPLOT),'   perturbative order: ', num2str(korder-1)]);
  %elseif (iopPLOT==5)
  %  title([num2str(latsize),'  bc',bc,'  kappa: 0.',kappa,'   operator: \psi \psi','   perturbative order: ', num2str(korder-1)]);
  %end
%  xlim([0,160]);
  %ylim([2.5*10^(-4),17*10^(-4)]);%iop=1:4
  %ylim([1*10^(-3),3.5*10^(-2)]);%for ioplot=5
  %ylim([4.5*10^(-4),1.8*10^(-3)]);%for ioplot=6:9 interpret
  %ylim([6,8]);
%   ylim([2.9*10^(-5),1.15*10^(-4)]);% iop=1:4
  %ylim([8.395*10^(-5),2.7*10^(-4)]);% iop=5:9 


%  filename = [num2str(latsize),'-',bc,'-',kappa,'-',num2str(iopPLOT),'.jpg'];
 % saveas(fig,filename);
  %hold off
%elseif(1==0)
 % hold on
  %plot(xmax,psERRAVG(:,6+1,iopPLOT),'-*b','LineWidth',1);
  %legend('0.1570','0.1560','0.1550')
  %hold off
%elseif(1==0)
 % hold on
  %plot(xmax,psERRAVG(:,4+1,iopPLOT),'-*r','LineWidth',1);
  %plot(xeigH,esgpsERRAVG(:,4+1,iopPLOT),'-dr','LineWidth',1);
  %plot(xmax,psERRAVG(:,6+1,iopPLOT),'-*g','LineWidth',1);
  %plot(xeigH,esgpsERRAVG(:,6+1,iopPLOT),'-dg','LineWidth',1);
  %hold off
%elseif(1==0)
 % hold on

 % if (korder - 1 == 1)
 %   plot(xmax,psERRAVG(:,iopPLOT),'r','LineWidth',3);
  %elseif (korder - 1 == 2)
   % plot(xmax,psERRAVG(:,iopPLOT),'b','LineWidth',3);
  %elseif (korder - 1 == 3)
   % plot(xmax,psERRAVG(:,iopPLOT),'g','LineWidth',3);
 % elseif (korder - 1 == 4)
  %  plot(xmax,psERRAVG(:,iopPLOT),'k','LineWidth',3);
  %end

  %title('Perturbative Error Results for Various Orders of \kappa within J_4');
  %legend('1^{st} Order \kappa','2^{nd} Order \kappa','3^{rd} Order \kappa','4^{th} Order \kappa');
  %hold off
end


C=((A/B)^2 -1)*100

clear val ii iop ieig xmax xeig xeigH nmax
