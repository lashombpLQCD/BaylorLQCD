%cd /home/qcd2/matlabtests/hermitian/trace/;
%matlabpool close force local
%matlabpool
%parpool(2)

%Here is my edit

  latsize = 12121216;   % lattice dimensions (string)
  k = 0.15700;       % kappa value
  ntrials = 1;     % trials
 % nrhs = 2; 
  nrhs = 10; %Commented out for testing -PL      % noises per trial % Chris changed from 200 -> 3
  bc = '1111';      % boundary conditions
  numconfigs = 10;
 % polys = [7]; 
  %polys = [50, 100, 200]; % Commented out for testing -PL
  polys = [50, 100, 200]; 
  % dimension of quark matrix
  if (latsize == 4444)
  	rank = 3072;
  elseif(latsize == 8888)
    rank = 49152;
  elseif(latsize == 12121216)
    rank = 331776
  end
  
  % no longer a need to declare nsub TW 4/9/2020
  %nsub = 199;         % maximum subtraction level in PE, deg-1 TW

  hvol = rank/12;   % hyper volume of lattice

% Chris -- loops over configUSED with cc
for polyrun = 1:length(polys)
	
	p = polys(polyrun)
	for runnum = 1:numconfigs  
  		configUSED = runnum;   % config used

		%%%%%% DISPLAY RUN INFORMATION %%%%%%%
  		disp(['Lattice Size:        ',num2str(latsize)]);
  		disp(['Configuration:       ',num2str(configUSED)]);
  		disp(['Boundary Conditions: ',bc]);
  		disp(['Rank:                ',num2str(rank)]);
  		disp(['k:                   ',num2str(k)]);
  		disp('');
  		disp(['Number of Trials:    ',num2str(ntrials)]);
  		disp(['Number of Z4 noise:  ',num2str(nrhs)]);
  		disp(['Seed:                ',num2str(1256+(configUSED-1)*rank*nrhs*ntrials+1)]);
  		disp(['Working Directory:   ',pwd]);
  		disp(['Storage Directory:   ','/data/lashombp/allsub_10noises/',num2str(latsize),'/bc',bc,'/',num2str(k*10000)]);

  		% file extension for config; {'01','02',...}
  		if (configUSED < 10)
    		ext = ['0',num2str(configUSED)];
  		else
    		ext = num2str(configUSED);
  		end
		
		% making directory BS 5/29.014
		mkdir(['/data/lashombp/allsub_10noises/',num2str(latsize),'/bc',bc,'/',num2str(k*10000)])

  		% move to proper director for final storage
  		% rand('seed',1256+(configUSED-1)*rank*nrhs*ntrials+1); % Chris
  
		% Chris -- creates random noise for each cc
		% Chris -- end noise making
		load(['/home/lashombp/revampedMatlab2/z2--', num2str(ext), '.mat']);

		% Chris -- edit file openers/loaders/ for corect directory and make them change with cc
		load(['/home/lashombp/revampedMatlab2/h', num2str(ext), '.mat']);
		
		latsize
		completetrace
  		cd(['/data/lashombp/allsub_10noises/',num2str(latsize),'/bc',bc,'/',num2str(k*10000)])

  		% removed saves to variables that no longer exist due to removal of esgpolyERR, esgpolyAVG, espsERR and espsAVG TW 4/9/2020
  		%save(['fig',ext,'-',num2str(latsize),'.mat'],'nsgAVG','nsgERR','esAVG','esERR','esgAVG',...
  		%                                             'esgERR','esgpolyAVG','esgpolyERR','esgpolylatestERR','esgpolylatestAVG','polyERR','polyAVG','psAVG','psERR',...
  		%                                             'espsAVG','espsERR','esgpsAVG','esgpsERR','neig','neigH','JaveHFES','JaveHFESPS','invgxi','invgxipoly')
		%
 		% save(['fig',ext,'-',num2str(latsize),'.mat'],'nsgAVG','nsgERR','esAVG','esERR','esgAVG',...
  		%                                             'esgERR','esgpolylatestERR','esgpolylatestAVG','polyERR','polyAVG','psAVG','psERR',...
   		%                                            'espsAVG','esgpsERR','neig','neigH','JaveHFES','JaveHFESPS','invgxi','invgxipoly')

 		save(['fig',ext,'-',num2str(latsize),'_deg',num2str(p),'_',num2str(nrhs),'noises.mat'],'nsgAVG','nsgERR','esAVG','esERR','esgAVG',...
                                               'esgERR','esgpolylatestERR','esgpolylatestAVG','polyERR','polyAVG','psAVG','psERR',...
                                               'esgpsERR','neig','neigH','JaveHFES','JaveHFESPS','invgxi','invgxipoly')

  		save(['vecmodes',ext,'-',num2str(latsize),'.mat'],'z','gx','gd','gvr','M','vr','vl','d');
		cd(['/home/lashombp/revampedMatlab2'])
		%delete(gcp('nocreate'))

	end %runconfig

end %polyrun


%save(['heigvals.mat'],'heig');

%save(['evectors.mat'],'veig');

%matlabpool close force local
%parpool close force local
%p = gcp;
%delete(gcp('nocreate'))
