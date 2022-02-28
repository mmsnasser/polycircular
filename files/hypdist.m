function [dis,zi,w0i] = hypdist (et,etp,n,alpha,z,w0)
% This function computes the hyperbolic distance dis between a point w0 (it
% is possible for w0 to be a vector of points too) and a row vector of 
% points z, in a simply connected domain G where:
% et, etp:  the parametrization of the boundary of G and its derivative 
% n: the number of discretization points
% alpha: a given point in G 
rhoD = @(x,y)(2*asinh(abs(x-y)./sqrt((1-abs(x).^2).*(1-abs(y).^2))));
A      = et-alpha;
gam    =-log(abs(et-alpha));
[mun,h]= fbie(et,etp,A,gam,n,5,[],1e-14,200);
fet    =(gam+h+i*mun)./A;
c      = exp(-mean(h(1:n)));
Phi    = c.*(et-alpha).*exp(gam+h+i.*mun);
%
z = z(:).'; w0 = w0(:).';
zi     = fcau(et,etp,Phi,z);
w0i    = fcau(et,etp,Phi,w0);
zi = zi(:); w0i = w0i(:);
dis=[];
for k=1:length(w0i)
    dis=[dis rhoD(zi,w0i(k))];
end
% 
end