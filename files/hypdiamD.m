function hd = hypdiamD(zet)
% This function computes the hyperbolic diameter of a comapct set E with
% respect to the unit disk D where zet is a row vector of points on the
% boundary of E
rhoD = @(x,y)(2*asinh(abs(x-y)./sqrt((1-abs(x).^2).*(1-abs(y).^2))));
z = zet(:); w = zet(:);
dis=[];
for k=1:length(w)
    dis=[dis rhoD(z,w(k))];
end
%
hd  = max(max(dis));
end