% makefig3b.m
clear
addpath ../bie; addpath ../fmm; addpath ../files
% 
n =   2^13
p0   = [-8   -4i    8    8i   -8];
d0   = [ 1    1];
[et0,et0p] = cirarcp3pt(p0,d0,n/2);
p1   = [-2   -4+i  -6   -4+2i -2];
d1   = [ 1   -1];
[et1,et1p] = cirarcp3pt(p1,d1,n/2);
p2   = [ 6    4+i   2    2+2i  6];
d2   = [ 1   -1];
[et2,et2p] = cirarcp3pt(p2,d2,n/2);
et  = [et0  ; et1  ; et2 ];  
etp = [et0p ; et1p ; et2p];
alphav = [-4+1.5i ; 4+2i ];
alpha = 0; m =2;
%
format long
cap = capac(et,etp,alphav,m,alpha)
abs(cap-10.95831564100293)
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
%
hold on;box on
plot(real(et0),imag(et0),'-k','LineWidth',1.5);
plot(real(et1),imag(et1),'-b','LineWidth',1.5);
plot(real(et2),imag(et2),'-b','LineWidth',1.5);
%
plot(real(p0(1:2:end)),imag(p0(1:2:end)),'sr','Markerfacecolor','r');
plot(real(p0(2:2:end)),imag(p0(2:2:end)),'dr','Markerfacecolor','r');
%
plot(real(p1(1:2:end)),imag(p1(1:2:end)),'sr','Markerfacecolor','r');
plot(real(p1(2:2:end)),imag(p1(2:2:end)),'dr','Markerfacecolor','r');
%
plot(real(p2(1:2:end)),imag(p2(1:2:end)),'sr','Markerfacecolor','r');
plot(real(p2(2:2:end)),imag(p2(2:2:end)),'dr','Markerfacecolor','r');
%
% plot(real(alphav),imag(alphav),'pr','Markerfacecolor','r');
% plot(real(alpha),imag(alpha),'pm','Markerfacecolor','r');
%
grid on; grid('minor')
set(gca, 'XMinorTick','on'); set(gca, 'YMinorTick','on')
ax=gca; ax.GridAlpha=0.5; ax.MinorGridAlpha=0.5;
set(gca,'FontSize',22)
set(gca,'LooseInset',get(gca,'TightInset'))
axis equal
axis([-9  9 -5  10])
xticks([-9:3:9])
print  -depsc sec-3ex-cirp-m2
%%