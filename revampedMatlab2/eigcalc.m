 function [heig, veig] = eigcalc(vl, vr, d, A)

%program to calculate the eigenvalues of the hermitian wilson matrix using the
%eigenmodes of the non hermitian wilson matrix

%normalization of left and right vectors
%only normalization of lefts is necessary in this fashion
%this preserves the orthogonality of the rights to themselves
for i = 1:2:size(vr,2)-1
    j = i+1;
    z = vl(:,i)'*vr(:,i);
    s(i) = z/sqrt((conj(z)*z));
    s(j) = conj(s(i));
    vl(:,i) = (s(i)*vl(:,i))/sqrt(sqrt(conj(z)*z));
    vl(:,j) = (s(j)*vl(:,j))/sqrt(sqrt(conj(z)*z));
    vr(:,i) = vr(:,i)/sqrt(sqrt(conj(z)*z));
    vr(:,j) = vr(:,j)/sqrt(sqrt(conj(z)*z));
end
    

tmp1 = zeros(1,size(vr,2));

for i = 1:2:size(vr,2)-1
    tmp1(i) = norm(vr(:,i));
end

tmp2 = zeros(1,size(vr,2));
for i = 2:2:size(vr,2)
    tmp2(i) = norm(vr(:,i));
end

%remove extraneous zero entries
for i = size(tmp1,2):-2:1
    tmp1(i) = [];
end

%remove extraneous zero entries
for i = size(tmp2,2)-1:-2:1
    tmp2(i) = [];
end

%tmp1 and tmp2 now hold the magnitudes of the V1R and V2R vectors
%respectively

%get mag of eval
for i = 1:size(d,1)
    mag(i) = d(i)'*d(i);
    mag(i) = sqrt(mag(i));
end

%get magnitude of s
for i = 1:size(s,2)
    smag(i) = s(i)'*s(i);
    smag(i) = sqrt(smag(i));
end

i = 1;
for i = 1:size(tmp1,2)
    tmp(i) = tmp1(i)*mag(i*2);
% end
% 
% 
% for i = 1:size(tmp2,2)
    tmptmp(i) = tmp2(i)*tmp(i);
% end
% 
% i = 1;
% for i = 1:size(tmptmp,2)
    tmp(i) = smag(i*2)*tmptmp(i); %removing s TW 9/25/17
% end
% 
% for i = 1:size(tmptmp,2)
    tmpsq(i) = tmp(i)'*tmp(i);
end


% V1R* dot V2R
for i = 1:2:size(vr,2)-1
%     for j = i+1:2:size(vr,2)
        j = i+1;
        dotprod(i) = vr(:,i)'*vr(:,j); %changed conjvr to vr TW 9/18/17 see DR. W's derivation
%     end
end

%remove zeros from column vector
for i = size(dotprod,2)-1:-2:1
    dotprod(i) = [];
end

%multiply the dotproduct by the eval
j = 2;
for i = 1:size(dotprod,2)
    tmpprod(i) = d(j).'*dotprod(i);
    j = j+2;
end

%multiply the above term by s
j = 1;
for i = 1:size(tmpprod,2)
    tmp(i) = s(j).'*tmpprod(i);
    j = j+2;
end

for i = 1:size(tmp,2)
    itmpprod(i) = imag(tmp(i))'*imag(tmp(i));
end

for i = 1:size(tmpsq,2)
    tmpdiff(i) = tmpsq(i)-itmpprod(i);
% end
% 
% for i = 1:size(tmpdiff,2)
    secterm(i) = sqrt(tmpdiff(i));
end

%for the first term of the eval
for i = 1:2:size(vr,2)-1
%     for j = i+1:2:size(vr,2)
        j = i+1;
        dotprod(i) = vr(:,i)'*vr(:,j);
%     end
end

for i = size(dotprod,2)-1:-2:1
    dotprod(i) = [];
end

for i = 1:size(dotprod,2)
    first(i) = d(i*2).'*dotprod(i);
end

j = 1;
for i = 1:size(first,2)
    first(i) = s(j).'*first(i);
    j = j+2;
end

for i = 1:size(first,2)
    first(i) = real(first(i));
end

for i = 1:size(first,2)
    Eplus(i) = first(i) + secterm(i);
end

for i = 1:size(first,2)
    Eminus(i) = first(i) - secterm(i); %changed to minus 9/19/17 TW to match Dr. W's derivation
end

heig = zeros(size(d,1),1);
j = 1;
for i = 1:size(Eplus,2)
    heig(j) = Eplus(i);
    j = j+2;
end

j = 2;
for i = 1:size(Eminus,2)
    heig(j) = Eminus(i);
    j = j+2;
end

[vals,indx] = sort(abs(heig));
heig = heig(indx);

%indice = 1:size(vr,2);
%ind = 1:size(gd);

%scatter(indice,heig)
%hold on;
%scatter(ind,gd)

clear i j  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% NEW SECTION %%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%generate the eigenvectors%%%

%calculate the phase of c2
for k = 1:size(tmp,2)
    imagtmp(k) = 1i*imag(tmp(k));
end

j = 1;
for i = 1:size(secterm,2)
    phaseplus(i) =(s(j)*d(j))*(secterm(i) + imagtmp(i))/(mag(j)*mag(j)*tmp1(i));
    phaseminus(i) =(s(j)*d(j))*(-1*secterm(i) + imagtmp(i))/(mag(j)*mag(j)*tmp1(i));
    j = j+2;
end

%calculate the mag of c2
for i = 1:size(secterm,2)
    taupp1(i) = (first(i)/secterm(i)) + 1;
    taump1(i) = ((-1*first(i))/secterm(i)) + 1;
    taupp1(i) = sqrt((1/2)*taupp1(i));
    taump1(i) = sqrt((1/2)*taump1(i));
end

for i = 1:size(phaseplus,2)
    c2p(i) = taupp1(i)*phaseplus(i);
    c2m(i) = taump1(i)*phaseminus(i);
end

%calculate c1

for i = 1:size(taupp1,2)
    c1p(i) = tmp1(i)*taupp1(i);
    c1m(i) = tmp1(i)*taump1(i);
end

%assemble the eigenvectors
veig = zeros(size(vl,1),size(vl,2));
j = 1;
k = 2;
for i = 1:size(c1p,2)
    veig(:,j) = c1p(i)*vl(:,j) + c2p(i)*vl(:,k);
    j = j+2;
    k = k+2;
end

j = 2;
k = 1;
for i = 1:size(c1m,2)
    veig(:,j) = c1m(i)*vl(:,k) + c2m(i)*vl(:,j);
    j = j+2;
    k = k+2;
end

veig = veig(:,indx);

[Mg5] = Mgamma5(A,size(A,1),size(A,2));
for i = 1:size(veig,2)
    resnorm(i) = norm(Mg5*veig(:,i)-heig(i)*veig(:,i));
end


clear i j k c1m c1p c2m c2p taupp1 taump1 tmp tmp1 tmp2 tmpsq tmpdiff imagtmp 
clear mag s smag tmpprod dotprod itmpprod Eminus Eplus first secterm tmptmp    
    
