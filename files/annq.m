function [q,cap,Phi] = annq (et,etp,n,zz,z2,type)
% annq.m
% Nasser, June 10, 2019
% This function computes the inner radius q for the conformal mapping
% w=Phi(z) from doubly connected domain G onto the annulus q<|w|<1 and
% the capacity of G, cap(G)=2pi/log(1/q),  where:
% et, etp:  the parametrization of the boundary of G and its derivative 
% n: the number of discretization points
% zz: zz=alpha is a given point in $G$ for bounded G; 
%     and zz=z1 is a given point in the interior of the curve that will be
%     maped onto the unit circle for unbounded G
% z2: a given point interior to the curve that will be maped onto the 
% circle |w|=q for both cases of bounded and unbounded G.
% type='b' for bounded G; and type='u' for unbounded G
%
if type=='b' 
    alpha = zz; A = et-alpha; gam = -log(abs((et-z2)./(alpha-z2)));
elseif type=='u' 
    z1 = zz; A = ones(size(et));  gam = -log(abs((et-z2)./(et-z1)));
end
[mun,h] =  fbie(et,etp,A,gam,n,5,[],1e-14,100);
q     =  exp(mean(h(n+1:2*n))-mean(h(1:n)));
cap   =  2*pi/log(1/q);
if type=='b' 
    Phi = exp(-mean(h(1:n))).*((et-z2)./(alpha-z2)).*exp(gam+h+i.*mun);
elseif type=='u' 
    Phi = exp(-mean(h(1:n))).*((et-z2)./(et-z1)).*exp(gam+h+i.*mun);
end
end