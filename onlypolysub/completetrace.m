
% --- clear data holders ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear M gM z gx x gd gvr d vr vl invxi xh 
clear nsg_raw esg_raw ps_raw ns_raw es_raw esps_raw esgps_raw
clear nsAVG nsgAVG esAVG esgAVG psAVG espsAVG esgpsAVG
clear nsERR nsgERR esAVG esgERR psERR espsERR esgpsERR
clear neig neigH

% --- prep initial info ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kB = k*B;
I = speye(rank);
M = I - kB;
M = sparse(M);
clear B
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
makemodes; % Chris
% load eigen.mat; % Chris
gdinv = inv(diag(gd));
dinv = inv(diag(d));

% --- GMRES INFO ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rtolDR=1e-20;
rtolPj=1e-8;
cyclemax=300;
mDR=80; kDR=40;
%p = 200 %degree of pert and poly expansions, now set in evscript.m 
% --- Solve 1st rhs (results are thrown away) ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z = z2(:,1);
tic;
[gx,gvk,gdk] = gmresdr5(M,z,mDR,kDR,rtolDR,cyclemax);
toc;

% Generate 1/xi for ESPS method, no longer needed, relic of the past 4/9/2020 TW
%*************************************************
%%invxi = zeros(neig,neig,nsub+1);
%Mp = vr;
%invxi(:,:,1) = diag(diag(vl' * Mp));
%for ksub = 2:nsub+1
%  Mp = vr + (kB * Mp);
%  invxi(:,:,ksub) = diag(diag(vl' * Mp));
%end
%******************************************


% Generate 1/gxi for gESPS method aka HFPS 
%invgxi = zeros(neigH,neigH,nsub+1); -PL
invgxi = sparse(neigH,neigH); % -PL

ggvr = gvr; %these are the eigenvectors of M' coming from makemodes
Mp = ggvr;
% BS 11/30/2014 this line replaced by line below it invgxi(:,:,1) = diag(diag(gvr' * Mp));

%invgxi(:,:,1) = diag(diag(gvr' * gammaX(Mp,5)));
invgxi(:,:) = diag(diag(gvr' * gammaX(Mp,5)));


for ksub = 2:p  % Removing loop over ksub -PL
  Mp = ggvr + (kB * Mp);
  %  BS 11/30/2014 this line replaced by line below it invgxi(:,:,ksub) = diag(diag(gvr' * Mp));
%invgxi(:,:,ksub) = diag(diag(gvr' * gammaX(Mp,5))); -PL  %the last term of this is all that is needed, when ksub = nsub+1!!!! -TW 4/9/2020
invgxi(:,:) = diag(diag(gvr' * gammaX(Mp,5))); % -PL 
end
clear Mp ggvr


%***first take at rewriting the above section ^^^^ TW 4/9/2020, I think this is ok
%Generate 1/gxi for gESPS method aka HFPS
%invgxi = sparse(neigH,neigH);
%Mp = gvr;
%for i = 1:deg %deg to be defined in evscript, deg = nsub +1, still need for loop to form the pert expansion TW 4/9/2020
%Mp = gvr + kB * Mp; %this will form the whole pert expansion 
%end
%invgxi = diag(diag(gvr' * gammaX(Mp,5))); %multiplies by gamma 5 on the right and gvr' on the left. Nested diags form a diagonal matrix
%clear Mp ggvr
%***********************************************

% invgxi % Chris 

%For invgxipoly
% load z2.mat; 
% load h01.mat;
%I=speye(49152);
%k = 0.1570;
%M = I-(k*B);
v = zeros(rank,1);
v(:,1)= z2(:,1);

%p = 200; % BS i think p is equivalent to ksub changed here TW
%n = 49152;
%for j=1:p
 %      v(:,j+1) = M * v(:,j);
 %  end
%lsmat = v(:,2:p+1)'*v(:,2:p+1);
%cls = v(:,2:p+1)'*v(:,1);
%  dls = lsmat\cls;
%  A = dls

%[~,rv,~] = NewtonBasisPP(M,v(:,1),p,[]);
[~,~,~,~,~,~,~,th] = gmresdrEIGritz(M,v(:,1),p,1,1e-2,1);
[rv] = ModLejaComplex(th);
%invgxipoly = zeros(neigH,neigH,nsub+1); % -PL
invgxipoly = sparse(neigH,neigH); % -PL
ggvr = gvr;
[phi] = newpoly(M,ggvr,p,rv);
%Mp = A(1,1)*ggvr;
%temp=gvr;




%invgxipoly(:,:,1) = diag(diag(gvr' * gammaX(phi(:,:,1),5))); -PL
%invgxipoly(:,:) = diag(diag(gvr' * gammaX(phi(:,:,1),5))); % -PL %getting rid of the last index TW 5/13/20 
invgxipoly(:,:) = diag(diag(gvr' * gammaX(phi(:,:),5)));
%for ksub = 2:nsub+1 Removing for loop over ksub -PL  
%temp=M*temp;
%Mp=Mp+A(ksub,1)*temp;
%invgxipoly(:,:,ksub) = diag(diag(gvr' * gammaX(Mp,5)));

%invgxipoly(:,:,ksub) = diag(diag(gvr' * gammaX(phi(:,:,ksub),5))); % -PL 

%end -PL





% --- Loop for various trials ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iruns=1:ntrials
  disp([' ---------- Run: ',num2str(iruns)]);

  % --- noise ---

%%%  z = z4noise(rank,nrhs);
%  load z2.mat; % Chris -- we now generate z2 in evscript
  %znoise(:,(iruns-1)*nrhs+1:iruns*nrhs) = z2;
  z = z2(:,1:nrhs);
  % --- Prep ES data with noise
  invgd_gvrdag_z = gdinv * (gvr' * z);
  invd_vldag_z = dinv * (vl' * z);

  % --- Prep PS data with noise




  %Minvpert_z = zeros(rank,nrhs,nsub+1); % -PL
  Minvpert_z = sparse(rank,nrhs);  % -PL
  %Minvpert_z(:,:,1) = z; %this doesn't work due to N-dimensional indexing clash with sparse call TW 5/12/20
  Minvpert_z(:,:) = z;
  for ksub=1:p
    %Minvpert_z(:,:,ksub+1) = z + (kB * Minvpert_z(:,:,ksub)); %same reason as above TW
     Minvpert_z(:,:) = z + (kB * Minvpert_z(:,:));
  end




% Now comes the part for polynomial method BS  4/3/2015
% Value of A was initially calculated elsewhere
%A = [   6.5614 - 0.0151i
 %     -21.0752 + 0.0696i
  %     40.3994 - 0.1562i
   %   -48.9310 + 0.2083i
    %   37.7279 - 0.1729i
     % -17.9445 + 0.0880i
      %  4.7977 - 0.0252i
       %-0.5516 + 0.0031i];
% Determine coefficients
% Chris -- comments out section mini-version
%load h01.mat
%I=speye(49152);
%k = 0.1570;
%M = I-(k*B);
% end Chris comment storm

% Chris -- changed from 
% v = zeros(rank,1);
% v(:,1)=z2(:,1);

%% NEED THIS FOR HIGH ORDER SUBTRACTION TEST - TW 12/1/19
A = zeros(rank,1);
A = z2(:,1);
v(:,1)= A;

%p = 200; % BS i think p is equivalent to ksub changed here, put p = degree above and commenting out all other declarations TW
%n = rank;
%for j=1:p
  %     v(:,j+1) = M * v(:,j);
 %  end
%lsmat = v(:,2:p+1)'*v(:,2:p+1);  
%cls = v(:,2:p+1)'*v(:,1);
  %dls = lsmat\cls;
  %A = dls

%[~,rv,~] = NewtonBasisPP(M,v(:,1),p,[]);
[~,~,~,~,~,~,~,th] = gmresdrEIGritz(M,v(:,1),p,1,1e-2,1);
[rv] = ModLejaComplex(th);
[phi] = newpoly(M,z,p,rv);

%Minvper_z = zeros(rank,nrhs,nsub+1); % -PL 
Minvper_z = sparse(rank,nrhs); % -PL

  % Minvper_z(:,:,1) = A(1,:)*z;
%Minvper_z(:,:,1) = phi(:,:,1);

%tryn = zeros(rank,nrhs,nsub+1); -PL
%tryn(:,:,1) = z;



  %for ksub=1:nsub+1 % Eliminating the for loop over ksub  -PL  %was 1, now 2 TW

    %tryn(:,:,ksub+1) = (M * tryn(:,:,ksub));

    %Minvper_z(:,:,ksub+1) =  Minvper_z(:,:,ksub)+A(ksub+1,:)* tryn(:,:,ksub+1);
    
	%Minvper_z(:,:,ksub) =  phi(:,:,ksub); % -PL
	%Minvper_z(:,:) =  phi(:,:,1); % -PL, this was needed in previous calc, but not here TW 5/12/20
  Minvper_z(:,:) = phi;
  %end % -PL
% END NEED THIS SECTION -TW 12/1/19




 % --- Prep ESPS  and gESPS data with noise
  %invxi_vldag_z = zeros(neig,nrhs,nsub+1); %for ESPS, relic of the past TW 4/9/2020
  

  %invgxi_gvrdag_z = zeros(neigH,nrhs,nsub+1); % -PL
  invgxi_gvrdag_z = sparse(neigH,nrhs); % -PL
  %invgxipoly_gvrdag_z = zeros(neigH,nrhs,nsub+1); % -PL
  invgxipoly_gvrdag_z = sparse(neigH,nrhs); % -PL
 

  %for ksub=1:nsub+1 % Removing for loop over ksub -PL
   % invxi_vldag_z(:,:,ksub) = invxi(:,:,ksub) * (vl' * z); %for ESPS relic of the past, TW 4/9/2020 
    %invgxi_gvrdag_z(:,:,ksub) = invgxi(:,:,ksub) * (gvr' * z); % -PL
    invgxi_gvrdag_z(:,:) = invgxi(:,:) * (gvr' * z); % -PL
    %invgxipoly_gvrdag_z(:,:,ksub) = invgxipoly(:,:,ksub) * (gvr' * z); % -PL
    invgxipoly_gvrdag_z(:,:) = invgxipoly(:,:) * (gvr' * z); % -PL
  %end


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

  nsg_raw   = zeros(nrhs,5); % NEEDED FOR HIGHER ORDER POLY SUB -TW 12/1/19
  es_raw    = zeros(nrhs,neig,5);
  esg_raw   = zeros(nrhs,neigH,5);
  
  %ps_raw    = zeros(nrhs,nsub+1,5); -PL
  ps_raw    = zeros(nrhs,5); % -PL
  %poly_raw  = zeros(nrhs,nsub+1,5); -PL
  poly_raw  = zeros(nrhs,5); % -PL
  
  %esps_raw  = zeros(nrhs,neig,nsub+1,5); %relic of the past TW 4/9/2020
  
  %esgps_raw = zeros(nrhs,neigH,nsub+1,5); -PL
  esgps_raw = zeros(nrhs,neigH,5); % -PL
  
  JaveHFES  = zeros(nrhs,neigH,5);

  %JaveHFESPS = zeros(nrhs,neigH,nsub+1,5); -PL
  JaveHFESPS = zeros(nrhs,neigH,5); % -PL

  %esgpoly_raw= zeros(nrhs,neigH,nsub+1,5); %incorrect way for HFPOLY TW 4/9/2020
  
  %esgpolylatest_raw= zeros(nrhs,neigH,nsub+1,5); -PL  %correct way for HFPOLY TW 4/9/2020
  esgpolylatest_raw= zeros(nrhs,neigH,5); % -PL  %correct way for HFPOLY TW 4/9/2020
  
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

    % --- ESg (gamma5) aka HFES ---
    xp = zeros(rank,1);
    for ieig=1:neigH
      xp = xp + invgd_gvrdag_z(ieig,rhs)*gvr(:,ieig);
      for iop = 1:5
        if (iop == 5) % scalar
          esg_raw(rhs,ieig,iop) = nsg_raw(rhs,iop) - z(:,rhs)'*gammaX(xp,5);
          JaveHFES(rhs,ieig,iop) = z(:,rhs)'*gammaX(xp,5);
        else % gamma_mu
          esg_raw(rhs,ieig,iop) = nsg_raw(rhs,iop) - z(:,rhs)'*gammaX(gammaX(xp,5),iop);
          JaveHFES(rhs,ieig,iop) =  z(:,rhs)'*gammaX(gammaX(xp,5),iop);
% esg_raw(rhs,ieig,iop) = nsg_raw(rhs,iop) - z(:,rhs)'*gammaX(gammaX(xp,iop),5);
 % BD 11/14/2014 interchanged position of iop and 5
       
       end
      end
    end

    % --- ES ---
    xp = zeros(rank,1);
    for ieig=1:neig
      xp = xp + invd_vldag_z(ieig,rhs)*vr(:,ieig);
      for iop = 1:5
        if (iop == 5) % scalar
          es_raw(rhs,ieig,iop) =  nsg_raw(rhs,iop) - z(:,rhs)'*xp;
        else % gamma_mu
          es_raw(rhs,ieig,iop) = nsg_raw(rhs,iop) - z(:,rhs)'*gammaX(xp,iop);
        end
      end
    end



    % --- PS ---
    %for ksub=1:nsub+1   Removing ksub For loop -PL
      for iop=1:5
        if (iop == 5) % scalar
          %ps_raw(rhs,ksub,iop) = nsg_raw(rhs,iop) - z(:,rhs)'*Minvpert_z(:,rhs,ksub); -PL
          ps_raw(rhs,iop) = nsg_raw(rhs,iop) - z(:,rhs)'*Minvpert_z(:,rhs); 
% is something missing in it?      
  else  % gamma_mu
          %ps_raw(rhs,ksub,iop) = nsg_raw(rhs,iop) - z(:,rhs)'*gammaX(Minvpert_z(:,rhs,ksub),iop); -PL
          ps_raw(rhs,iop) = nsg_raw(rhs,iop) - z(:,rhs)'*gammaX(Minvpert_z(:,rhs),iop);
% --- clear data holders ---
% --- clear data holders ---
% --- clear data holders ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
      end
    %end


% Now is the polynomial part BS 4/3/2015
  %for ksub=1:nsub+1    Removing ksub For loop -PL 
      for iop=1:5
        if (iop == 5) % scalar
          %poly_raw(rhs,ksub,iop) = nsg_raw(rhs,iop) - z(:,rhs)'*Minvper_z(:,rhs,ksub); -PL
          poly_raw(rhs,iop) = nsg_raw(rhs,iop) - z(:,rhs)'*Minvper_z(:,rhs);
% is something missing in it?
  else  % gamma_mu
          %poly_raw(rhs,ksub,iop) = nsg_raw(rhs,iop) - z(:,rhs)'*gammaX(Minvper_z(:,rhs,ksub),iop); -PL
          poly_raw(rhs,iop) = nsg_raw(rhs,iop) - z(:,rhs)'*gammaX(Minvper_z(:,rhs),iop);
% --- clear data holders ---
% --- clear data holders ---
% --- clear data holders ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
      end
  %end








    % --- ESPS ---, this needs to be commented out, relic of the past -TW 4/9/2020
%**************************************************************************
   % xp = zeros(rank,nsub+1);
   % for ieig = 1:neig
   %   for ksub = 1:nsub+1
   %     xp(:,ksub) = xp(:,ksub) + invxi_vldag_z(ieig,rhs,ksub)*vr(:,ieig);
   %     for iop=1:5
   %       if (iop == 5) % scalar
   %         esps_raw(rhs,ieig,ksub,iop) = es_raw(rhs,ieig,iop) + ps_raw(rhs,ksub,iop) - nsg_raw(rhs,iop) + (z(:,rhs)' * xp(:,ksub));
   %       else
   %         esps_raw(rhs,ieig,ksub,iop) = es_raw(rhs,ieig,iop) + ps_raw(rhs,ksub,iop) - nsg_raw(rhs,iop) + (z(:,rhs)' * gammaX(xp(:,ksub),iop));
         % should ksub and iop be interchanged?
   %        end
   %     end
   %   end
   % end
%****************************************************************************

    % --- ESgPS aka HFPS ---
    %xp = zeros(rank,nsub+1);
    xp = zeros(rank,1);
    for ieig = 1:neigH
      %for ksub = 1:nsub+1  - Removed for loop over ksub -PL 
        %xp(:,ksub) = xp(:,ksub) + invgxi_gvrdag_z(ieig,rhs,ksub)*gvr(:,ieig); -PL
       % xp(:) = xp(:) + invgxi_gvrdag_z(ieig,rhs)*gvr(:,ieig); %this line threw an error on the vector addition for some reason TW 5/13/20
        xp = xp + invgxi_gvrdag_z(ieig,rhs)*gvr(:,ieig); %replaced to be like ES subtraction TW 5/13/20
		for iop=1:5
          if (iop == 5) % scalar
            %esgps_raw(rhs,ieig,ksub,iop) = esg_raw(rhs,ieig,iop) + ps_raw(rhs,ksub,iop) - nsg_raw(rhs,iop) + (z(:,rhs)' * gammaX(xp(:,ksub),5)); -PL
            %JaveHFESPS(rhs,ieig,ksub,iop) = (z(:,rhs)' * gammaX(xp(:,ksub),5)); -PL
            esgps_raw(rhs,ieig,iop) = esg_raw(rhs,ieig,iop) + ps_raw(rhs,iop) - nsg_raw(rhs,iop) + (z(:,rhs)' * gammaX(xp(:),5));
            JaveHFESPS(rhs,ieig,iop) = (z(:,rhs)' * gammaX(xp(:),5));
		  else
            %esgps_raw(rhs,ieig,ksub,iop) = esg_raw(rhs,ieig,iop) + ps_raw(rhs,ksub,iop) - nsg_raw(rhs,iop) + (z(:,rhs)' * gammaX(gammaX(xp(:,ksub),5),iop)); -PL
            %JaveHFESPS(rhs,ieig,ksub,iop) =  (z(:,rhs)' * gammaX(gammaX(xp(:,ksub),5),iop)); -PL
            esgps_raw(rhs,ieig,iop) = esg_raw(rhs,ieig,iop) + ps_raw(rhs,iop) - nsg_raw(rhs,iop) + (z(:,rhs)' * gammaX(gammaX(xp(:),5),iop));
            JaveHFESPS(rhs,ieig,iop) =  (z(:,rhs)' * gammaX(gammaX(xp(:),5),iop));
		  end
        end
      %end
    end


% This is the part where i worked on HFPOLY BS 5/6/2015, this needs to be commented out, the error is calculated incorrectly. It uses the values coming from pert expansion -TW 4/9/2020
%******************************************************************
   % xp = zeros(rank,nsub+1);
   % for ieig = 1:neigH
   %   for ksub = 1:nsub+1
   %     xp(:,ksub) = xp(:,ksub) + invgxi_gvrdag_z(ieig,rhs,ksub)*gvr(:,ieig); %this step here is wrong. invgxi_gvrdag_z is from pert expansion TW 4/9/2020 
   %     for iop=1:5
   %       if (iop == 5) % scalar
   %         esgpoly_raw(rhs,ieig,ksub,iop) = esg_raw(rhs,ieig,iop) + poly_raw(rhs,ksub,iop) - nsg_raw(rhs,iop) + (z(:,rhs)' * gammaX(xp(:,ksub),5));
   %       else
   %         esgpoly_raw(rhs,ieig,ksub,iop) = esg_raw(rhs,ieig,iop) + poly_raw(rhs,ksub,iop) - nsg_raw(rhs,iop) + (z(:,rhs)' * gammaX(gammaX(xp(:,ksub),5),iop));
   %       end
   %     end
   %   end
   % end
%*****************************************************************



% This is the part where i worked on HFPOLYLATEST BS 5/6/2015, this is the correct version TW 4/9/2020
    % -- HFPOLY aka ESgPOLY -- 
    %xp = zeros(rank,nsub+1); -PL
    xp = zeros(rank,1);
    for ieig = 1:neigH
      %for ksub = 1:nsub+1 - Removing for loop over ksub -PL 
        %xp(:,ksub) = xp(:,ksub) + invgxipoly_gvrdag_z(ieig,rhs,ksub)*gvr(:,ieig); -PL
        %xp(:) = xp(:) + invgxipoly_gvrdag_z(ieig,rhs)*gvr(:,ieig); %weird thing with array dims -TW 5/13/120
        xp = xp + invgxipoly_gvrdag_z(ieig,rhs)*gvr(:,ieig);
        for iop=1:5
          if (iop == 5) % scalar
            %esgpolylatest_raw(rhs,ieig,ksub,iop) = esg_raw(rhs,ieig,iop) + poly_raw(rhs,ksub,iop) - nsg_raw(rhs,iop) + (z(:,rhs)' * gammaX(xp(:,ksub),5)); -PL
            esgpolylatest_raw(rhs,ieig,iop) = esg_raw(rhs,ieig,iop) + poly_raw(rhs,iop) - nsg_raw(rhs,iop) + (z(:,rhs)' * gammaX(xp(:),5));
          else
            %esgpolylatest_raw(rhs,ieig,ksub,iop) = esg_raw(rhs,ieig,iop) + poly_raw(rhs,ksub,iop) - nsg_raw(rhs,iop) + (z(:,rhs)' * gammaX(gammaX(xp(:,ksub),5),iop)); -PL
            esgpolylatest_raw(rhs,ieig,iop) = esg_raw(rhs,ieig,iop) + poly_raw(rhs,iop) - nsg_raw(rhs,iop) + (z(:,rhs)' * gammaX(gammaX(xp(:),5),iop));
          end
        end
      %end
    end














    toc;
  end


 %start correcting for the real and imaginary results BS 12/8/2015


% Previous dimensions when using ksub for loops -PL 

%{
            for iop = 1:5
                if (iop ==5)
                   esgpolylatest_raw(:,:,:,iop) = real(esgpolylatest_raw(:,:,:,iop)); %correct HFPOLY TW
                   %esgpoly_raw(:,:,:,iop) = real(esgpoly_raw(:,:,:,iop)); %incorrect HFPOLY TW
                   esgps_raw(:,:,:,iop)   = real(esgps_raw(:,:,:,iop));
                   %esps_raw(:,:,:,iop)    = real(esps_raw(:,:,:,iop)); %relic of the past TW
                   poly_raw(:,:,iop)      = real(poly_raw(:,:,iop)); 
                   ps_raw(:,:,iop)        = real(ps_raw(:,:,iop));
                   es_raw(:,:,iop)        = real(es_raw(:,:,iop));
                   esg_raw(:,:,iop)       = real(esg_raw(:,:,iop));
                   nsg_raw(:,iop)         = real(nsg_raw(:,iop));
                   JaveHFES(:,:,iop)    = real(JaveHFES(:,:,iop));
                   JaveHFESPS(:,:,:,iop)    = real(JaveHFESPS(:,:,:,iop));        
               else
                   esgpolylatest_raw(:,:,:,iop) = imag(esgpolylatest_raw(:,:,:,iop)); %correct HFPOLY TW
                   %esgpoly_raw(:,:,:,iop) = imag(esgpoly_raw(:,:,:,iop)); %incorrect HFPOLY TW
                   esgps_raw(:,:,:,iop)   = imag(esgps_raw(:,:,:,iop));
                   %esps_raw(:,:,:,iop)    = imag(esps_raw(:,:,:,iop)); %relic of the past TW
                   poly_raw(:,:,iop)      = imag(poly_raw(:,:,iop)); 
                   ps_raw(:,:,iop)        = imag(ps_raw(:,:,iop));
                   es_raw(:,:,iop)        = imag(es_raw(:,:,iop));
                   esg_raw(:,:,iop)       = imag(esg_raw(:,:,iop));
                   nsg_raw(:,iop)         = imag(nsg_raw(:,iop));
                   JaveHFES(:,:,iop)      = imag(JaveHFES(:,:,iop));
                   JaveHFESPS(:,:,:,iop)  = imag(JaveHFESPS(:,:,:,iop));
              end
          end
%}


            for iop = 1:5
                if (iop ==5)
                   esgpolylatest_raw(:,:,iop) = real(esgpolylatest_raw(:,:,iop)); %correct HFPOLY TW
                   %esgpoly_raw(:,:,:,iop)    = real(esgpoly_raw(:,:,:,iop)); %incorrect HFPOLY TW
                   esgps_raw(:,:,iop)         = real(esgps_raw(:,:,iop));
                   %esps_raw(:,:,:,iop)       = real(esps_raw(:,:,:,iop)); %relic of the past TW
                   poly_raw(:,iop)            = real(poly_raw(:,iop)); 
                   ps_raw(:,iop)              = real(ps_raw(:,iop));
                   es_raw(:,:,iop)            = real(es_raw(:,:,iop));
                   esg_raw(:,:,iop)           = real(esg_raw(:,:,iop));
                   nsg_raw(:,iop)             = real(nsg_raw(:,iop));
                   JaveHFES(:,:,iop)          = real(JaveHFES(:,:,iop));
                   JaveHFESPS(:,:,iop)        = real(JaveHFESPS(:,:,iop));        
               else
                   esgpolylatest_raw(:,:,iop) = imag(esgpolylatest_raw(:,:,iop)); %correct HFPOLY TW
                   %esgpoly_raw(:,:,:,iop)    = imag(esgpoly_raw(:,:,:,iop)); %incorrect HFPOLY TW
                   esgps_raw(:,:,iop)         = imag(esgps_raw(:,:,iop));
                   %esps_raw(:,:,:,iop)       = imag(esps_raw(:,:,:,iop)); %relic of the past TW
                   poly_raw(:,iop)            = imag(poly_raw(:,iop)); 
                   ps_raw(:,iop)              = imag(ps_raw(:,iop));
                   es_raw(:,:,iop)            = imag(es_raw(:,:,iop));
                   esg_raw(:,:,iop)           = imag(esg_raw(:,:,iop));
                   nsg_raw(:,iop)             = imag(nsg_raw(:,iop));
                   JaveHFES(:,:,iop)          = imag(JaveHFES(:,:,iop));
                   JaveHFESPS(:,:,iop)        = imag(JaveHFESPS(:,:,iop));        
              end
          end



  % Avg over rhs's
  for iop=1:5
   % nsg_raw_AVG(iruns,iop) = mean(nsg_raw(:,iop))/hvol;
 nsgAVG(iop) = mean(nsg_raw(:,iop))/hvol;
 nsgERR(iop) = std(nsg_raw(:,iop))/(hvol*sqrt(nrhs));

    for ieig=1:neigH
     % esg_raw_AVG(iruns,ieig,iop) = mean(esg_raw(:,ieig,iop))/hvol;
      esgAVG(ieig,iop) = mean(esg_raw(:,ieig,iop))/hvol;
  esgERR(ieig,iop) = std(esg_raw(:,ieig,iop))/(hvol*sqrt(nrhs));
 
 
 %for ksub=1:nsub+1  Removing the for loop over ksub -PL 
      %  esgps_raw_AVG(iruns,ieig,ksub,iop) = mean(esgps_raw(:,ieig,ksub,iop)) / hvol;

% esgpsAVG(ieig,ksub,iop) = mean(esgps_raw(:,ieig,ksub,iop)) / hvol; -PL
% esgpsERR(ieig,ksub,iop) = std(esgps_raw(:,ieig,ksub,iop)) / (hvol*sqrt(nrhs)); -PL
% JaveHFESPS(ieig,ksub,iop)=std(JaveHFESPS(:,ieig,ksub,iop)) / (hvol*sqrt(nrhs)); -PL
 esgpsAVG(ieig,iop) = mean(esgps_raw(:,ieig,iop)) / hvol;
 esgpsERR(ieig,iop) = std(esgps_raw(:,ieig,iop)) / (hvol*sqrt(nrhs));
 JaveHFESPS(ieig,iop) = std(JaveHFESPS(:,ieig,iop)) / (hvol*sqrt(nrhs));
% This part is added for Hermitian forced polynomial subtraction BS 5/6/2015
    


%esgpolylatestAVG(ieig,ksub,iop) = mean(esgpolylatest_raw(:,ieig,ksub,iop)) / hvol; -PL
esgpolylatestAVG(ieig,iop) = mean(esgpolylatest_raw(:,ieig,iop)) / hvol;

%***********%esgpolyAVG and esgpolyERR calculated incorrectly, esgpolylatest are the ones to keep, TW 4/9/2020
%esgpolyAVG(ieig,ksub,iop) = mean(esgpoly_raw(:,ieig,ksub,iop)) / hvol;



 %esgpolylatestERR(ieig,ksub,iop) = std(esgpolylatest_raw(:,ieig,ksub,iop)) / (hvol*sqrt(nrhs)); -PL  
 esgpolylatestERR(ieig,iop) = std(esgpolylatest_raw(:,ieig,iop)) / (hvol*sqrt(nrhs));  

% esgpolyERR(ieig,ksub,iop) = std(esgpoly_raw(:,ieig,ksub,iop)) / (hvol*sqrt(nrhs));  

 %end %ksub


    end %ieig

    for ieig=1:neig
 %     es_raw_AVG(iruns,ieig,iop) = mean(es_raw(:,ieig,iop))/hvol;
 esAVG(ieig,iop) = mean(es_raw(:,ieig,iop))/hvol;
 esERR(ieig,iop) = std(es_raw(:,ieig,iop))/(hvol*sqrt(nrhs));

%ESPS is a relic of the past, commenting out as it is not necessary TW 4/9/2020
      %for ksub=1:nsub+1
       % esps_raw_AVG(iruns,ieig,ksub,iop) = mean(esps_raw(:,ieig,ksub,iop)) / hvol;
 %espsAVG(ieig,ksub,iop) = mean(esps_raw(:,ieig,ksub,iop)) / hvol;
 %espsERR(ieig,ksub,iop) = std(esps_raw(:,ieig,ksub,iop)) /(hvol*sqrt(nrhs));
     
 %end %ksub

    end %ieig

    %for ksub=1:nsub+1  Commented out for loop over ksub -PL 
  %    ps_raw_AVG(iruns,ksub,iop) = mean(ps_raw(:,ksub,iop))/hvol;
       
       %psAVG(ksub,iop) = mean(ps_raw(:,ksub,iop))/hvol; -PL
       %psERR(ksub,iop) = std(ps_raw(:,ksub,iop))/(hvol*sqrt(nrhs)); -PL
       psAVG(iop) = mean(ps_raw(:,iop))/hvol;
       psERR(iop) = std(ps_raw(:,iop))/(hvol*sqrt(nrhs));
    %end

% Added part for poly BS 4/3/2015 NEED THIS -TW 12/1/19
    %for ksub=1:nsub+1  Commented out for loop over ksub -PL 
  %    ps_raw_AVG(iruns,ksub,iop) = mean(ps_raw(:,ksub,iop))/hvol;
       
      %polyAVG(ksub,iop) = mean(poly_raw(:,ksub,iop))/hvol; -PL 
      %polyERR(ksub,iop) = std(poly_raw(:,ksub,iop))/(hvol*sqrt(nrhs)); -PL 
       polyAVG(iop) = mean(poly_raw(:,iop))/hvol;
       polyERR(iop) = std(poly_raw(:,iop))/(hvol*sqrt(nrhs));
    %end



  end %IOP
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

clear znoise gxsol xsol;
clear iruns rhs ksub xp xpp ieig iop;
clear Minvpert_z cyclemax eig dinv dk ;
clear gdinv gdk gvk invd_vrdag_z invgd_gvrdag_z;

%clear znoise gxsol xsol;
%clear iruns rhs ksub xp xpp ieig iop;
%clear Minvpert_z cyclemax eig dinv dk es_raw es_raw_AVG;
%clear esps_raw_AVG esgps_raw_AVG;
%clear esg_raw esg_raw_AVG esps_raw esgps_raw;
%clear gdinv gdk gvk invd_vrdag_z invgd_gvrdag_z;

clear kB kDR mDR ;
%clear kB kDR mDR ns_raw ns_raw_AVG ps_raw ps_raw_AVG;
