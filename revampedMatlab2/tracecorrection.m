numconfigs = 10; 
latsize = 8888;   % lattice dimensions (string)
k = 0.15700; %kappa
bc = '1111';
p = 10; %degree of polynomial

for runnum = 1:numconfigs
configUSED = runnum;
%configUSED = 1;
  if (configUSED < 10)
    ext = ['0',num2str(configUSED)];
  else
    ext = num2str(configUSED);
  end


load(['/home/lashombp/revampedMatlab2/h', num2str(ext), '.mat']); %load hopping matrix
load(['/home/lashombp/revampedMatlab2/z2--', num2str(ext), '.mat']); %load rhs from subtraction
z = z2(:,1); %this needs to be the same rhs used to generate the ritz values in the subtraction part
n = size(B);

M = speye(n)-k*B; %form the Wilson Dirac operator
clear B

unitvecs = speye(n); %this is all the unit vectors needed to form the trace

nc = 3; ns = 4; %color, spin indices
nx = 8; ny = 8; nz = 8; nt = 8; %space, time indices
npx = 2; npy = 2; npz = 2; npt = 1; %numprocs in each direction
numprocs = npx*npy*npz*npt; %total number of processors in each direction, must match requested number in run_trace_correction.sh

ix = nx/npx; iy = ny/npy; iz = nz/npz; it = nt/npt; %ix iy iz and it must be an integer

[~,~,~,~,~,~,~,th] = gmresdrEIGritz(M,z,p,1,1e-2,1); %calculate ritz values
[rv] = ModLejaComplex(th); %order them according to modified leja

tic()
parpool('local',numprocs) %start parallel pool of workers
spmd(numprocs)
  e_i = unitvecs(:,(labindex-1)*ix*iy*iz*it*nc*ns+1:labindex*ix*iy*iz*it*nc*ns); %making the indices explicit, this splits up the unit vectors and sends them to each processor
  scalarterms = tracecorrection_scalar(M,e_i,p,rv); %all of the terms from each processor
  localterms = tracecorrection_local(M,e_i,p,rv);
end
toc()

scalartermscell = scalarterms(:); %sends the terms to the client from the workers
localtermscell = localterms(:);

cd(['/data/lashombp/revampedMatlab2/',num2str(latsize),'/bc',bc,'/',num2str(k*10000)])
save(['scalarcorrectionterms.mat'],'scalartermscell');
save(['localcorrectionterms.mat'],'localtermscell');

delete(gcp('nocreate'))

%gather and sum all the terms
scalarcorrection = 0;
for i = 1:numprocs
  scalarcorrection = scalarcorrection + cell2mat(scalartermscell(i));
end
scalarcorrection
save(['scalarcorrection_config',num2str(ext),'.mat'],'scalarcorrection')

localcorrection = zeros(ns,1);
for i = 1:numprocs
  tmp = cell2mat(localtermscell(i));
  for mu = 1:ns
    localcorrection(mu,1) = localcorrection(mu,1) + tmp(mu,1);
  end
end
localcorrection
save(['localcorrection_config',num2str(ext),'.mat'],'localcorrection')

end %loop over configs
