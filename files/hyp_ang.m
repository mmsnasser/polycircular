function beta2 = hyp_ang(b1,b2,b3)
% This function compute the angles beta2 between the hyperbolic sides 
% b2b1 and b2b3.
T        =  @(z,a)((z-a)./(1-conj(a).*z));
omg      =   abs(carg(T(b3,b2)/T(b1,b2)));
beta2    =   Arg(T(b3,b2),angle(T(b1,b2)))-angle(T(b1,b2));
end