
% --- clear data holders ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
clear M gM z gx x gd gvr d vr vl invxi xh 
clear nsg_raw esg_raw ps_raw ns_raw es_raw esps_raw esgps_raw
clear nsAVG nsgAVG esAVG esgAVG psAVG espsAVG esgpsAVG
clear nsERR nsgERR esAVG esgERR psERR espsERR esgpsERR
clear neig neigH
%}


clear M z x 
clear nsg_raw poly_raw 
clear nsgAVG polyAVG
clear nsgERR polyERR



% --- prep initial info ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kB = k*B;
I = speye(rank);
M = I - kB;
M = sparse(M);
%clear B
latsize
% --- Holder for noise and solution vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ntrials = 1;
% nrhs = 200; % Chris -- defined in evscript
 % BS 5/23/014 znoise = zeros(rank,ntrials * nrhs);
 % BS 5/23/014 gxsol  = zeros(rank,ntrials * nrhs);
 % BS 5/23/014 xsol   = zeros(rank,ntrials * nrhs);

% --- EIG INFO --%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % BS 5/23/014 makemodes;
%makemodes; % Chris
% load eigen.mat; % Chris
%gdinv = inv(diag(gd));
%dinv = inv(diag(d));

% --- GMRES INFO ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rtolDR=1e-20;
rtolPj=1e-8;
cyclemax=300;
mDR=80; kDR=40;

% --- Solve 1st rhs (results are thrown away) ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z = z2(:,1);
tic;
[gx,gvk,gdk] = gmresdr5(M,z,mDR,kDR,rtolDR,cyclemax);
toc;



%v = zeros(rank,1);
%v(:,1)= z2(:,1);

%[~,~,~,~,~,~,~,th] = gmresdrEIGritz(M,v(:,1),p,1,1e-2,1);
%[rv] = ModLejaComplex(th);




% --- Loop for various trials ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iruns=1:ntrials
  	disp([' ---------- Run: ',num2str(iruns)]);


  	z = z2(:,1:nrhs); 

  	A = zeros(rank,1); 
  	A = z2(:,1); 
   	v(:,1) = A; 
  
  	[~,~,~,~,~,~,~,th] = gmresdrEIGritz(M,v(:,1),p,1,1e-2,1);
  	[rv] = ModLejaComplex(th);
  	[phi] = newpoly(M,z,p,rv);

  	Minvper_z = sparse(rank,nrhs); % -PL

  	Minvper_z(:,:) = phi;


  	% --- Solution Vectors ---
  	gx = zeros(rank,nrhs);
  	parfor rhs=1:nrhs
    	tic;
    	disp(['rhs: ',num2str(rhs)]);
    	gx(:,rhs) = gmresproj5(M,z(:,rhs),mDR,kDR,gvk,gdk,rtolPj,cyclemax);
    	toc;
  	end
  	gxsol(:,(iruns-1)*nrhs+1:iruns*nrhs) = gx;


  	% --- Run RHS ---

	nsg_raw = zeros(nrhs,5); 
  	poly_raw  = zeros(nrhs,5); % -PL
  
  	for rhs=1:nrhs
    	tic; 
    	disp([' ----- Rhs: ',num2str(rhs)]);
		
		% --- NS ---
		for iop=1:5
      		if (iop == 5)  % scalar
        		nsg_raw(rhs,iop) = z(:,rhs)'*gammaX(gx(:,rhs),5);
      		else  % gamma_mu
        		nsg_raw(rhs,iop) = z(:,rhs)'*gammaX(gammaX(gx(:,rhs),5),iop);
 				% BD 11/14/2014 interchanged position of iop and 5
      		end
    	end

		% --- POLY --- 
      	for iop=1:5
        	if (iop == 5) % scalar
          		poly_raw(rhs,iop) = nsg_raw(rhs,iop) - z(:,rhs)'*Minvper_z(:,rhs);
	   			%Javepoly(rhs,iop) = z(:,rhs)'*phi(:,rhs); Don't think I need -PL
        	else  % gamma_mu
          		poly_raw(rhs,iop) = nsg_raw(rhs,iop) - z(:,rhs)'*gammaX(Minvper_z(:,rhs),iop);
	  			%Javepoly(rhs,iop) = z(:,rhs)'*gammaX(phi(:,rhs),iop); Don't think I need -PL
				% --- clear data holders ---
				% --- clear data holders ---
				% --- clear data holders ---
				%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        	end
      	end

      	for iop = 1:5
         	if (iop ==5)
            	poly_raw(:,iop)         = real(poly_raw(:,iop));
	    		%Javepoly(:,iop)         = real(Javepoly(:,iop)); Don't think I need -PL
			else
            	poly_raw(:,iop)         = imag(poly_raw(:,iop));
	    		%Javepoly(:,iop) 	   = imag(Javepoly(:,iop)); Don't think I need -PL
         	end
      	end


  		for iop=1:5
       		nsgAVG(iop) = mean(nsg_raw(:,iop))/hvol; 
			nsgERR(iop) = std(nsg_raw(:,iop))/(hvol*sqrt(nrhs)); 
		
			polyAVG(iop) = mean(poly_raw(:,iop))/hvol;
			polyERR(iop) = std(poly_raw(:,iop))/(hvol*sqrt(nrhs));
			%Jave(iop) = mean(Javepoly(:,iop))/hvol; Don't think I need -PL
		end %IOP

 	end %nrhs
end %iruns

% Average over all trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for iop=1:5
%  if (iop==5) % scalar (use real)
%    nsgAVG(:,iop) = ones(neig,1) * mean(nsg_raw_AVG(:,iop));
%    nsgERR(:,iop) = ones(neig,1) * std(real(nsg_raw_AVG(:,iop)))/sqrt(ntrials);

 %   for ieig=1:neigH
  %    esgAVG(ieig,iop) = mean(esg_raw_AVG(:,ieig,iop));
   %   esgERR(ieig,iop) = std(real(esg_raw_AVG(:,ieig,iop)))/sqrt(ntrials);
    %  for ksub=1:nsub+1
     %   esgpsAVG(ieig,ksub,iop) = mean(esgpnns_raw_AVG(:,ieig,ksub,iop));
      %  esgpsERR(ieig,ksub,iop) = std(real(esgps_raw_AVG(:,ieig,ksub,iop)))/sqrt(ntrials);
      %end
   % end

   % for ieig=1:neig
   %   esAVG(ieig,iop) = mean(es_raw_AVG(:,ieig,iop));
    %  esERR(ieig,iop) = std(real(es_raw_AVG(:,ieig,iop)))/sqrt(ntrials);
     % for ksub=1:nsub+1
      %  espsAVG(ieig,ksub,iop) = mean(esps_raw_AVG(:,ieig,ksub,iop));
       % espsERR(ieig,ksub,iop) = std(real(esps_raw_AVG(:,ieig,ksub,iop)))/sqrt(ntrials);
      %end
   % end

   % for ksub=1:nsub+1
    %  psAVG(ksub,iop) = mean(ps_raw_AVG(:,ksub,iop));
     % psERR(ksub,iop) = std(real(ps_raw_AVG(:,ksub,iop)))/sqrt(ntrials);
   % end
 % else % gamma_mu (use imag)
  %  nsgAVG(:,iop) = ones(neig,1) * mean(nsg_raw_AVG(:,iop));
   % nsgERR(:,iop) = ones(neig,1) * std(imag(nsg_raw_AVG(:,iop)))/sqrt(ntrials);

   % for ieig=1:neigH
    %  esgAVG(ieig,iop) = mean(esg_raw_AVG(:,ieig,iop));
     % esgERR(ieig,iop) = std(imag(esg_raw_AVG(:,ieig,iop)))/sqrt(ntrials);
      %for ksub=1:nsub+1
       % esgpsAVG(ieig,ksub,iop) = mean(esgps_raw_AVG(:,ieig,ksub,iop));
        %esgpsERR(ieig,ksub,iop) = std(imag(esgps_raw_AVG(:,ieig,ksub,iop)))/sqrt(ntrials);
     % end
    %end

    %for ieig=1:neig
     % esAVG(ieig,iop) = mean(es_raw_AVG(:,ieig,iop));
      %esERR(ieig,iop) = std(imag(es_raw_AVG(:,ieig,iop)))/sqrt(ntrials);
      %for ksub=1:nsub+1
       % espsAVG(ieig,ksub,iop) = mean(esps_raw_AVG(:,ieig,ksub,iop));
        %espsERR(ieig,ksub,iop) = std(imag(esps_raw_AVG(:,ieig,ksub,iop)))/sqrt(ntrials);
      %end
    %end

   % for ksub=1:nsub+1
    %  psAVG(ksub,iop) = mean(ps_raw_AVG(:,ksub,iop));
    %  psERR(ksub,iop) = std(imag(ps_raw_AVG(:,ksub,iop)))/sqrt(ntrials);
    %end
 % end
%end
% Store all noise/solution vectors
z=z2(:,1:nrhs); gx=gxsol;


%{
clear znoise gxsol xsol;
clear iruns rhs ksub xp xpp ieig iop;
clear Minvpert_z cyclemax eig dinv dk ;
clear gdinv gdk gvk invd_vrdag_z invgd_gvrdag_z;
%}
%clear znoise gxsol xsol;
%clear iruns rhs ksub xp xpp ieig iop;
%clear Minvpert_z cyclemax eig dinv dk es_raw es_raw_AVG;
%clear esps_raw_AVG esgps_raw_AVG;
%clear esg_raw esg_raw_AVG esps_raw esgps_raw;
%clear gdinv gdk gvk invd_vrdag_z invgd_gvrdag_z;

clear kB kDR mDR ;
%clear kB kDR mDR ns_raw ns_raw_AVG ps_raw ps_raw_AVG;
