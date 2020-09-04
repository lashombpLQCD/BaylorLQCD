function y = applyppgamma5(inA,invec,n,deg,g,diracAlg,basis,rval,gam)

%===============================================================
%   A function to implement the polynomial preconditioner p(A) by 
%   applying it to a vector z and returning it (i.e. y = p(A)z).
%   Function not restricted to left or right preconditioning
%===============================================================
%   INPUT: 
%       
%       inA     = Matrix in p(A) and in Ax=b
%       invec   = Vector z to apply p(A) to
%       n       = Size of A
%       deg     = Polynomial p(A) is degree (deg-1)
%       g       = Coefficients of polynomial p(A)
%       diracAlg= if set to string "useDiracAlg" it will use Dirac Algebra
%                 gamma matrices in mvps. 
%       basis   = Basis used to construct p(A)
%       rval    = Modified Leja ordered Ritz values from Lanczos 
%                 for Newton Basis polynomial p(A) 
%       gam     = Variable used in construction of Newton Basis 
%                 polynomial p(A)
%
%----------------------------------------------------------------
%   OUTPUT:
%
%       y       = invec with polynomial preconditioner p(A) 
%                 applied to it. 
%---------------------------------------------------------------

y = zeros(n,1); 
w = zeros(n,deg);  

switch basis

    case 'POWER'
        w(:,1) = invec(:,1); 
        y(:,1) = g(1)*w(:,1);

        for i = 1:deg-1
            if strcmp(diracAlg,"useDiracAlg")
                w(:,i+1) = inA*gamma5(w(:,i),n,1);
            else
                w(:,i+1) = inA*w(:,i);
            end
            pGlobalmvps; 
            y(:,1) = y(:,1) + g(i+1)*w(:,i+1);
        end 
    %End of Power
    
    case 'NEWTON'
        w(:,1) = invec(:,1);
        y(:,1) = g(1)*w(:,1);
        
        j = 1; 
        while j < deg
            if strcmp(diracAlg,"useDiracAlg")
                t(:,1) = inA*gamma5(w(:,j),n,1);
            else
                t(:,1) = inA*w(:,j);
            end
            pGlobalmvps;
                
            %NEED TO FIX. IF HERMITIAN, SHOULD BE ZERO BUT STILL PICKING UP
            %A SMALL IMAGINARY PART
            if(imag(rval(j)) < 1e-8 ) %If Ritz values are real
                t(:,1) = t(:,1) - rval(j)*w(:,j); 
                w(:,j+1) = t(:,1)/gam(j+1); 
                y(:,1) = y(:,1) + g(j+1)*w(:,j+1); 
                j = j + 1; 
            else %If Ritz values are complex 
                t(:,1) = t(:,1) - real(rval(j))*w(:,j);
                if strcmp(diracAlg,"useDiracAlg")
                    s(:,1) = inA*gamma5(t(:,1),n,1);  
                else
                    s(:,1) = inA*t(:,1);
                end
                pGlobalmvps;
                s(:,1) = s(:,1) - real(rval(j))*t + (imag(rval(j)))^2*w(:,j); 
                w(:,j+1) = t(:,1)/gam(j+1); y(:,1) = y(:,1) + g(j+1)*w(:,j+1); 
                w(:,j+2) = s(:,1)/gam(j+2); y(:,1) = y(:,1) + g(j+2)*w(:,j+2); 
                j = j + 2; 
            end % End ifelse
        end % End While
    %End Newton
    
end %End Switch
end %End of function