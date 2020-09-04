function [g,rv,gam] = NewtonBasisPP(A,b,deg,diracAlg)
%==========================================================================
%   A Minres polynomial preconditioner using the Newton Basis for
%   Ritz values
%==========================================================================
%  INPUT:
%       A       = an nxn symmetric matrix
%       b       = a single right-hand side
%       deg     = polynomial is degree (deg-1)
%       diracAlg= if set to string "useDiracAlg" it will use Dirac Algebra
%                 gamma matrices in mvps. 
%--------------------------------------------------------------------------
%  OUTPUT:
%   y       = the initial residual vector for the polynomial preconditioned
%             system (i.e. r = b)
%   g       = coefficients for polynomial 
%--------------------------------------------------------------------------

n = size(A,1); 

r = b;
%transy(:,1) = r;
vp(:,1) = b/norm(b); 


%Use one cycle of Lanczos to obtain Harmonic Ritz values 

%Maybe I should just use gmresdrEIG5, instead, since it's small enough 
%that full orthogonalization won't be that expensive to matter

if(deg <= 1)
    %deg < 1 means that the polynomial is of degree 0 and won't need Ritz
    %values
    rv = 0; 
    gam = 1; 
else
    %rv = lanczosEIG5(A,b(:,1),deg-1,diracAlg);
    %rv = arnoldiEIG5(A,b(:,1),deg-1,diracAlg);
    [~,~,~,~,~,~,~,th] = gmresdrEIGritz(A,b,deg,1,1e-2,1); 
%     dabs = abs(ritz);
%     [~,ind] = sort(dabs);
%     ritz = ritz(ind);
    %Leja ordering the Ritz values
    %rv = ModLeja(ritz);
    [rv] = ModLejaComplex(th);

    %rv = leja(ritz);
    %Formulation of polynomial (assuming real Ritz values)
    gam = ones(deg-1);
end    
    
D = zeros(n,deg); 

j = 1; 

while j < deg 
    if strcmp(diracAlg,"useDiracAlg")
        t(:,1) = A*gamma5(vp(:,j),n,1); D(:,j) = t(:,1);  
    else
        t(:,1) = A*vp(:,j); D(:,j) = t(:,1);
    end
    %pGlobalmvps; 
    
    %NEED TO FIX THIS. SHOULD BE REAL, BUT HAVE VERY, VERY TINY IMAGINARY
    %PART
    if(imag(rv(j)) < 1e-8) %If Ritz values are real
        t(:,1) = t(:,1) - rv(j)*vp(:,j); gam(j+1) = norm(t); 
        vp(:,j+1) = t(:,1)/gam(j+1);  
        j = j + 1; 
    else %If Ritz values are complex 
        fprintf("in imag part\n")
        t(:,1) = t(:,1) - real(rv(j))*vp(:,j); gam(j+1) = norm(t); 
        if strcmp(diracAlg,"useDiracAlg")
            s(:,1) = A*gamma5(t(:,1),n,1); D(:,j+1) = s(:,1)/gam(j+1);  
        else
            s(:,1) = A*t(:,1); D(:,j+1) = s/gam(j+1);
        end
        %pGlobalmvps;
        s(:,1) = s(:,1) - real(rv(j))*t + (imag(rv(j)))^2*vp(:,j);
        vp(:,j+1) = t/gam(j+1); 
        gam(j+2) = norm(s); vp(:,j+2) = t(:,1)/gam(j+2); %POSSIBLY A TYPO %changed to t from s
        j = j + 2; 
   end 
        
end %End of while loop 

if strcmp(diracAlg,"useDiracAlg")
    D(:,deg) = A*gamma5(vp(:,deg),n,1);
else
    D(:,deg) = A*vp(:,deg); 
end 
%pGlobalmvps; 

%Solving D'Dg=D'v1 for g

for i = 1:deg
    for j = 1:deg
        beta = D(:,i)'*D(:,j);
   lsmat(i,j) = beta;
    end
end

for i=1:deg
    beta = D(:,i)'*vp(:,1);
    cls(:,i) = beta;
end

%Solving for g

g = lsmat\cls';

return 
