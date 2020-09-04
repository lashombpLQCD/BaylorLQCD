function [corrterm] = tracecorrection_local(M,e_i,p,rv)

%b = sparse(size(M,1),1);
%b(1) = 1;

%[~,~,~,~,~,~,~,th] = gmresdrEIGritz(M,b,p,1,1e-2,1);
%[rv] = ModLejaComplex(th);

corrterm = zeros(4,1);
for mu = 1:4
  for i = 1:size(e_i,2)
  %[~,~,~,~,~,~,~,th] = gmresdrEIGritz(M,e_i(:,i),p,1,1e-2,1);
  %[rv] = ModLejaComplex(th);
  [phi] = newpoly(M,e_i(:,i),p,rv);
  corrterm(mu,1) = corrterm(mu,1) + e_i(:,i)'*gammaX(phi,mu);
  end
end
