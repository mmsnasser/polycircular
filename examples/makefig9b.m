% makefig9b.m
% 02-06-2021
clear;
addpath ../bie; addpath ../fmm; addpath ../files; addpath ../pcm
%% 
nv   =  9*2.^[8:11];
%%
for kk=1:length(nv)
    n    =  nv(kk);
t    = [0:2*pi/n:2*pi-2*pi/n].';
% 
et0  =  exp(i.*t); 
et0p =  i.*exp(i.*t); 
%
ver1 = [ 3/9+1i/9  3/9-3i/9 -3/9-3i/9 -3/9+1i/9  4i/9 -1/9+1i/9  2/9+4i/9  1/9+1i/9  4/9+4i/9];
dir1 = [ 0         0         0        -1         1    -1         1        -1         1 ];

[ca3,~]=my3Pts(-3/9+1i/9,-2/9+3i/9,4i/9);   c3 = 4/9+1i/9;
[ca2,~]=my3Pts(-1/9+1i/9,3i/9,2/9+4i/9);    c2 = 6/9+1i/9;
[ca1,~]=my3Pts(1/9+1i/9,2/9+3i/9,4/9+4i/9); c1 = 8/9+1i/9;

slop2 = -(-1/9-real(ca2))/(1/9-imag(ca2));
tht2  =  (pi/2-atan(slop2))*180/pi;
slop1 = -(1/9-real(ca1))/(1/9-imag(ca1));
tht1  =  (pi/2-atan(slop1))*180/pi;
%%
cnt1 = [ inf  inf  inf  ca3  c3  ca2  c2  ca1  c1];
[et1,et1p] = plgsegcirarcp(ver1,cnt1,dir1,n/9);
%
deltv = [0;1];
et   = [et0 ;  et1  ];
etp  = [et0p ; et1p ];
alphav =  0.0i;
alpha  = -0.5i;
cap(kk,1)    = capm(et,etp,alphav,deltv,1,alpha);
%
end
%
cap
%%