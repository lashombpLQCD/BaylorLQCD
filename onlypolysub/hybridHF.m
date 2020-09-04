function [xh,thh,evectorh,hk,vk,grealh] = hybridHF(A,b,m,k,rtol,cyclim,fid)
%"Hybrid" inverter used to solve the hermitian forced case
%The linear equations of non hermitian case are calculated and used to solve for
%the linear equations of the hermitian case. The solution vector which minimizes 
% min ||c-Hd|| is found from gamma5*xm = x0 + Vm*d, and the algorithm proceeds as normal \
%from there

n = size(A,1);

if (nargin < 7)
  fid=1;
  if (nargin < 5)
    rtol = 1e-6;
    cyclim = 10;
    if (nargin < 3)
      m = 10;
      k = 5;
    end
  end
end


xnh = zeros(n,1);
cycle = 1;
j = 1;
jj = 1;

rninit = norm(b(:,1));
rn = rninit;
r = b;
vn = norm(r);
v(:,1) = r/vn;
c(1,1) = vn;
c(2:m+1,1) = zeros(m,1);

rnhinit = norm(b(:,1));
rnh = rnhinit;
rh = b;
vnh = norm(rh);
vh(:,1) = rh/vnh;
ch(1,1) = vnh;
ch(2:m+1,1) = zeros(m,1);

%while ( (rn/rninit > rtol) && (cycle <= cyclim) )
while ( (cycle <= cyclim) )
  while ( (j <= m) ) 
    wv = A*v(:,j);
    vnf = norm(wv);

    for i = 1:j
      h(i,j) = v(:,i)'*wv;
      wv = wv - h(i,j) * v(:,i);
    end
    vn = norm(wv);

    %--------reorthogonalization section-------%
    if( vn < 1.1*vnf )
      for i = 1:j
        dot = v(:,i)'*wv;
        wv = wv - dot * v(:,i);
        h(i,j) = h(i,j) + dot;
      end
      vn = norm(wv);
    end
    %--------------------------------------------------%

    h(j+1,j) = vn;
    v(:,j+1) = wv/h(j+1,j);

  % output data per cycle
  % disp(['iteration: ',num2str(j)]);

    j = j + 1;
  end
  
  j = 1;

  d(1:m,1) = h(1:m+1,1:m) \ c(1:m+1) ;
  srv(1:m+1,1) = c(1:m+1)-(h(1:m+1,1:m)*d(1:m,1));

  %----Set up and solve linear equations.-----%
  xnh(:,1) = xnh(:,1) + v(:,1:m)*d(1:m); 
  r = v(:,1:m+1)*srv;
  rn = norm(r);

  dspnhlineqns = ['Residual norm of non hermitian eqns = ', num2str(rn)];
  disp(dspnhlineqns);
  
 % if (cycle == 1000)
 %     a = xnh(:,1)
 % end
 
  hh = h(1:m,1:m);
  em = zeros(m,1);
  em(m) = 1;
  ff = hh' \ em;

  hh(:,m) = hh(:,m) + h(m+1,m)^2 * ff;
  hh=full(hh); [g,dd] = eig(hh,'nobalance'); 
% no balance needed!!!-WW; 03/25/2014!

  dd = diag(dd);
  dabs = abs(dd);
  [dabs,ind] = sort(dabs);
  th = dd(ind);
  gg = g(:,ind(1:k));
 
  for i=1:k
    rho(i) = gg(:,i)'*h(1:m,:)*gg(:,i);
    tv = h(1:m,:)*gg(:,i)-rho(i)*gg(:,i);
    tvn = norm(tv);
    rna(cycle,i) = sqrt( tvn*tvn+ abs(h(m+1,m))^2*abs(gg(m,i))^2 );
    tha(cycle,i) = th(i);
    rhoa(cycle,i) = rho(i);
  end

  [rnasort,ind] = sort(rna(cycle,:));
  rna(cycle,1:k)=rna(cycle,ind(1:k));
  gg = gg(:,ind(1:k));
  th = th(ind);

  greal = gg;

  greal(m+1,1:k) = zeros(1,k);


% Chris 

  beta = h(m+1,m);
  punty = zeros(m,1);
  punty = -beta*ff;
  greal(1:m,k+1) = punty(1:m,1);
  greal(m+1,k+1) = 1.0;

%  greal(:,k+1) = srv;

% end Chris


  [gon,rr] = qr(greal(:,1:k+1),0);
  hnew = gon'*h*gon(1:m,1:k);
  h(k+1,:) = zeros(1,m);

  j = 1;
  rtolev = 1e-11;
  while ( j<= k && rna(cycle,j) <= rtolev)
    hnew(j+1:k+1,j) = zeros(k-j+1,1);
    j = j + 1;
  end



  evector = v(:,1:m) * greal(1:m,1:k);


  h(1:k+1,1:k) = hnew;

  c(1:k+1,1) = gon(:,1:k+1)'*srv(:,1);
  c(k+2:m+1,1) = zeros(m-k,1);

  work = v*gon;
  v(:,1:k+1) = work;

  %section for just reorthog. one vector, v_{k+1}
  for i = 1:k
    dot = v(:,i)'*v(:,k+1) ;
    v(:,k+1) = v(:,k+1) - dot * v(:,i);
  end

  v(:,k+1) = v(:,k+1)/norm(v(:,k+1));

  % output data per cycle
  %fprintf(fid,'Cycle: %d  Rel Res Norm: %12.8g\n',cycle,rn/rninit);

  j = k + 1;
  cycle = cycle + 1;
end
  
  
  %************STARTING HF SECTION *******************%
cycle = 1;

rnhinit = norm(b(:,1));
rnh = rnhinit;
rh = b;
vnh = norm(rh);
vh(:,1) = rh/vnh;
ch(1,1) = vnh;
ch(2:m+1,1) = zeros(m,1);
xh = gamma5(xnh(:,1),n,1);

 while ( (cycle <= cyclim) )
  while ( (jj <= m) )
    wvh =A*gamma5(vh(:,jj),n,1); 
    vnfh = norm(wvh);

    for i = 1:jj
      H(i,jj) = vh(:,i)'*wvh; 
      wvh = wvh- H(i,jj) * vh(:,i);
    end
    vnh = norm(wvh);

    %--------reorthogonalization section-------%
    if( vnh < 1.1*vnfh )
      for i = 1:jj
        dot = vh(:,i)'*wvh;
        wvh = wvh - dot * vh(:,i);
        H(i,jj) = H(i,jj) + dot;
      end
      vnh = norm(wvh);
    end
    %--------------------------------------------------%

    H(jj+1,jj) = vnh;
    vh(:,jj+1) = wvh/H(jj+1,jj);

  % output data per cycle
  % disp(['iteration: ',num2str(j)]);

    jj = jj + 1;
 
  end

  %xnhg5(:,1) = gamma5(xnh(:,1),n,1);
  %for i = 1:m
   %   dotprod = vh(:,i)'*vh(:,i);
    %  if ( (dotprod > 0.9999) && (dotprod < 1.0001))
     %   dh(i,1) = vh(:,i)'*(xnhg5(:,1)-xh(:,1));%back solving for soln to least squares problem, TW
     % end
  %end

 % dspxnhg5 = ['first elements of xnhg5 = ', num2str(xnhg5(1)),' ', num2str(xnhg5(2)),' ', num2str(xnhg5(3))];
 % disp(dspxnhg5);
 % dspxh = ['first elements of xh = ', num2str(xh(1)),' ', num2str(xh(2)),' ', num2str(xh(3))];
 % disp(dspxh);

 
 
  dh(1:m,1) = H(1:m+1,1:m) \ ch(1:m+1);%testing for accuracy of dh, TW
  srvh(1:m+1,1) = ch(1:m+1)-(H(1:m+1,1:m)*dhtest(1:m,1));%change dtest back to dh, TW
  
  

  
  %----Set up and solve linear equations.-----%
%  xh(:,1) = gamma5(xnh(:,1),n,1); %added to get soln to lin eqns from non hermitian system TW
  xh(:,1) = xh(:,1) + vh(:,1:m)*dh(1:m);%for testing accuracy of xh, TW
 % xh(:,1) = gamma5(xnh(:,1),n,1); %added to get soln to lin eqns from non hermitian system TW
  %xh(:,1) = gamma5*xnh(:,1);
  %srvh(1:m+1,1) = ch(1:m+1)-(H(1:m+1,1:m)*dh(1:m,1));
  rh = vh(:,1:m+1)*srvh;
  rnh = norm(rh);

 % hnorm = norm(H);
 % vhnorm = norm(vh);
 % srvhnorm = norm(srvh);

 % dispresnorm = ['Residual norm of hermitian eqns = ',num2str(rnh)];
 % disphnorm = ['Norm of h matrix =',num2str(hnorm)];
 % dispvhnorm = ['Norm of vh = ',num2str(vhnorm)];
 % dispsrvnorm = ['Norm of srvh = ',num2str(srvhnorm)];

 % dspdh = ['first elements of dh = ', num2str(dh(1)),' ', num2str(dh(2)),' ', num2str(dh(3))];
 % dspch = ['first elements of ch = ', num2str(ch(1)),' ', num2str(ch(2)),' ', num2str(ch(3))];
 % dspsrvh = ['first elements of srvh = ', num2str(srvh(1)),' ', num2str(srvh(2)),' ', num2str(srvh(3))];
 % dspvh = ['first elements of vh = ', num2str(vh(1,1)),' ', num2str(vh(2,1)),' ', num2str(vh(3,1))];
  
  %testnorm = norm(dtest-dh);%accuracy testing TW
  %disp(testnorm);
 % linnorm = norm(xhtest-xh);%accuracy testing TW
  %disp(linnorm);

  %disptestnorm = ['Norm dif of dtest and dh = ',num2str(testnorm)];%print norm TW
 % displinnorm = ['Norm dif of xhtest and xh = ',num2str(linnorm)];%print norm TW

  %disp(disptestnorm);%print norm TW
 % disp(displinnorm);%print norm TW
  %disp(dispresnorm);
  %disp(disphnorm);
  %disp(dispvhnorm);
  %disp(dispsrvnorm);
  % disp(dspdh);
  % disp(dspsrvh);
  % disp(dspvh);
  % disp(dspch);

  %resnorm = norm(A*xh-b);
  dspresnorm = ['Residual norm of hermitian linear equations is = ', rnh];
  disp(dspresnorm);
  
  HH = H(1:m,1:m);
  em = zeros(m,1);
  em(m) = 1;
  ffh = HH' \ em;

  HH(:,m) = HH(:,m) + H(m+1,m)^2 * ffh;
  HH=full(HH); [gh,ddh] = eig(HH,'nobalance'); 
% no balance needed!!!-WW; 03/25/2014!

  ddh = diag(ddh);
  dabsh = abs(ddh);
  [dabsh,indh] = sort(dabsh);
  thh = dd(indh);
  ggh = gh(:,indh(1:k));
 
  for i=1:k
    rhoh(i) = ggh(:,i)'*H(1:m,:)*ggh(:,i);
    tvh = H(1:m,:)*ggh(:,i)-rhoh(i)*ggh(:,i);
    tvnh = norm(tvh);
    rnah(cycle,i) = sqrt( tvnh*tvnh+ abs(H(m+1,m))^2*abs(ggh(m,i))^2 );
    thah(cycle,i) = thh(i);
    rhoah(cycle,i) = rhoh(i);
  end
  


  [rnasorth,indh] = sort(rnah(cycle,:));
  rnah(cycle,1:k)=rnah(cycle,indh(1:k));
  ggh = ggh(:,indh(1:k));
  thh = thh(indh);

  grealh = ggh;

  grealh(m+1,1:k) = zeros(1,k);
  grealh(:,k+1) = srvh;

  [gonh,rrh] = qr(grealh(:,1:k+1),0);
  hnewh = gonh'*H*gonh(1:m,1:k);
  H(k+1,:) = zeros(1,m);

  jj = 1;
  rtolev = 1e-11;
  %rtolev = 1e-14;
  while (jj <= k && rnah(cycle,jj) <= rtolev)
    hnewh(jj+1:k+1,jj) = zeros(k-jj+1,1);
    jj = jj + 1;
  end



  evectorh = vh(:,1:m) * grealh(1:m,1:k);

  dspevec = ['Size of h evector is ',size(evectorh)];
  disp(dspevec);
  dspeval = ['Size of h eval is ',size(thh)];
  disp(dspeval);


  H(1:k+1,1:k) = hnewh;

  ch(1:k+1,1) = gonh(:,1:k+1)'*srvh(:,1);
  ch(k+2:m+1,1) = zeros(m-k,1);

  workh = vh*gonh;
  vh(:,1:k+1) = workh;

  %section for just reorthog. one vector, v_{k+1}
  for i = 1:k
    dot = vh(:,i)'*vh(:,k+1) ;
    vh(:,k+1) = vh(:,k+1) - dot * vh(:,i);
  end

  vh(:,k+1) = vh(:,k+1)/norm(vh(:,k+1));

  % output data per cycle
  %fprintf(fid,'Cycle: %d  Rel Res Norm: %12.8g\n',cycle,rn/rninit);

  jj = k + 1;
  cycle = cycle + 1;
end

hk = H(1:k+1,1:k);
vk = vh(:,1:k+1);

%xh;
%xhtest = gamma5*xnh;
% xhtest
 %dif = abs(dh-dtest)
 %ndif = norm(dif)
%dh;
%dtest;

%for i=1:k
% dispeig = ['value of eig ', i,' ', num2str(thh(i))];
 %disp(dispeig);
%end

if (rn/rninit < rtol) && (cycle-1 <= cyclim)
  fprintf(fid,'gmresdrEIG(%d,%d) converged in %d cycles with a relative residual of %12.8g\n', ...
          m,k,cycle-1,rn/rninit);
else
  fprintf(fid,'gmresdrEIG(%d,%d) stoped after %d cycles without converging to the desired tolerence of %12.8g\n', ...
          m,k,cycle-1,rn/rninit);
  fprintf(fid,'a relative residual of %12.8g has been reached\n',rn/rninit);
end

if (rnh/rnhinit < rtol) && (cycle-1 <= cyclim)
  fprintf(fid,'gmresdrEIG5(%d,%d) converged in %d cycles with a relative residual of %12.8g\n', ...
          m,k,cycle-1,rnh/rnhinit);
else
  fprintf(fid,'gmresdrEIG5(%d,%d) stoped after %d cycles without converging to the desired tolerence of %12.8g\n', ...
          m,k,cycle-1,rnh/rnhinit);
  fprintf(fid,'a relative residual of %12.8g has been reached\n',rnh/rnhinit);
end

return
