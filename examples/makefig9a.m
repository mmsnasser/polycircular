% makefig9a.m
% 02-06-2021
clear;
addpath ../bie; addpath ../fmm; addpath ../files; addpath ../pcm
%% 
n    =  9*2^7;
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
tht2  =  (pi/2-atan(slop2))*180/pi
slop1 = -(1/9-real(ca1))/(1/9-imag(ca1));
tht1  =  (pi/2-atan(slop1))*180/pi
%%
cnt1 = [ inf  inf  inf  ca3  c3  ca2  c2  ca1  c1];
[et1,et1p] = plgsegcirarcp(ver1,cnt1,dir1,n/9);
%
deltv = [0;1];
et   = [et0 ;  et1  ];
etp  = [et0p ; et1p ];
alphav =  0.0i;
alpha  = -0.5i;
cap    = capm(et,etp,alphav,deltv,1,alpha)
%
% [~,capa] = annq (et,etp,n,alpha,alphav,'b')
%%
figure;
hold on
box on
plot(real(et0),imag(et0),'-k','LineWidth',1.5);
plot(real(et1),imag(et1),'-k','LineWidth',1.5);
for k=4:2:8
    JJ=(k-1)*n/9+1:k*n/9;
    crv = et1(JJ);
    plot(real(crv),imag(crv),'-r','LineWidth',1.5);
end
md = [-2/9+3i/9 ; 3i/9 ; 2/9+3i/9];
% plot(real(md),imag(md),'dr','Markerfacecolor','r');
for k=5:2:9
    JJ=(k-1)*n/9+1:k*n/9;
    crv = et1(JJ);
    plot(real(crv),imag(crv),'-b','LineWidth',1.5);
end
cd = [4/9+1i/9; 6/9+1i/9; 8/9+1i/9];
% plot(real(cd),imag(cd),'pb','Markerfacecolor','b');
plot(real(ver1),imag(ver1),'sk','Markerfacecolor','k');
grid on; grid('minor')
set(gca, 'XMinorTick','on'); set(gca, 'YMinorTick','on')
ax=gca; ax.GridAlpha=0.5; ax.MinorGridAlpha=0.5;
set(gca,'FontSize',18)
axis equal
set(gca,'LooseInset',get(gca,'TightInset'))
axis([-1.1  1.1 -1.1  1.1])
print  -depsc BSimpson
print  -dpdf  BSimpson
drawnow
%%    