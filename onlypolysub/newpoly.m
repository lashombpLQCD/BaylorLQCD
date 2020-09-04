%old way looking at every sub level TW 4/9/2020
%function [phi] = newpoly(A,b,deg,rv);

%n = size(A,1);
%nb = size(b);

%poly = sparse(nb(1),nb(2));
%phi = zeros(nb(1),nb(2),deg);
%product = b;
%i = 1;

%while (i <= deg-1)
%    poly = poly + (1/rv(i))*product;
%    product = product - (1/rv(i))*A*product;
%    phi(:,:,i) = poly;
%    i = i + 1;
%end
%phi(:,:,i) = poly + (1/rv(i))*product;
%phi = poly;

%new way with reduced dimensionality of arrays TW 4/9/2020, leaving old way up for now
%**************************************************************
function [phi] = newpoly(A,b,deg,rv);

n = size(A,1);
nb = size(b);

poly = sparse(nb(1),nb(2));
phi = sparse(nb(1),nb(2));
product = b;
i = 1;

while (i <= deg-1)
    poly = poly + (1/rv(i))*product;
    product = product - (1/rv(i))*A*product;
    phi = poly;
    i = i + 1;
end
phi = poly + (1/rv(i))*product;
%**************************************************************

