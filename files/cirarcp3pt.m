function [et,etp] = cirarcp3pt(pt,dir,ns)
% cirarcp3.m
% Nasser, Dec 29, 2020
% This function compute the discretization of the parametrization of the
% circular arc polygon where, for k=1,2,...,m, the arc k passes through the
% three points: pt(2k-1), pt(2k), pt(2k+1), and
% and "ns" is the number of graded points on each arc of the polygon. 
% The number of graded points for the whole polygon is 
% "n=number of sides*ns". The graded mesh points are computed by the 
% function "deltw.m".
% 
if pt(end)==pt(1)
    ver =  pt(1:2:end-1);
    m   =  length(ver);
else
    pt(end+1)=pt(1);
    ver =  pt(1:2:end-1);
    m   =  length(ver);
end
%
for k=1:m
    [cent(k),rad(k)]=my3Pts(pt(2*k-1),pt(2*k),pt(2*k+1));
end
%
n =  m*ns;
t = (0:2*pi/n:2*pi-2*pi/n).';
%
[s,sp]     =   deltw(t,m,3);
%
ver(m+1) = ver(1);
for k=1:m
    a   = ver(k);   b = ver(k+1);  
    ang{k}    = angle([a-cent(k),b-cent(k)]);
    if dir(k)>0 & (ang{k}(2)<ang{k}(1))
        ang{k}(2)=ang{k}(2)+2*pi;
    elseif dir(k)<0 & (ang{k}(1)<ang{k}(2))
        ang{k}(1)=ang{k}(1)+2*pi;
    end
end
for k=1:m
    alp(k)=ang{k}(1);
    bet(k)=ang{k}(2);
end
for k=1:m
    thet = (m*(bet(k)-alp(k))/(2*pi)).*s(1+(k-1)*n/m:k*n/m)+k*alp(k)-(k-1)*bet(k);
    thetp= (m*(bet(k)-alp(k))/(2*pi)).*sp(1+(k-1)*n/m:k*n/m);
    et (1+(k-1)*n/m:k*n/m,1)= cent(k)+rad(k).*exp(i.*thet);
    etp(1+(k-1)*n/m:k*n/m,1)= rad(k).*exp(i.*thet).*(i.*thetp);
end
%
end