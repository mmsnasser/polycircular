function hd = hypdiam (et,etp,n,alpha,z)
% This function computes the hyperbolic diameter of a comapct set E 
% with respect to a simply connected domain G where:
% et, etp:  the parametrization of the boundary of G and its derivative 
% n: the number of discretization points
% alpha: a given point in G 
% z: is a row vector of points z in the boundary of E
hd  = max(max(hypdist(et,etp,n,alpha,z,z)));
end