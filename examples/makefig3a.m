% makefig3a.m
clear
addpath ../bie; addpath ../fmm; addpath ../files
% 
n =   2^13
p0   = [ 1+i  2i -1+i -2 -1-i -2i 1-i  2 1+i];
d0   = [ 1    1   1    1];
[et0,et0p] = cirarcp3pt(p0,d0,n/2);
p1   = [ 0.5+0.5i  0.4  0.5-0.5i -0.4i -0.5-0.5i -0.4 -0.5+0.5i  0.4i  0.5+0.5i];
d1   = [ 1         1    1         1];
[et1,et1p] = cirarcp3pt(p1,d1,n/2);
et  = [et0  ; et1  ];  
etp = [et0p ; et1p ];
alphav = [ 0 ];
alpha = i; m =1;
%
cap = capac(et,etp,alphav,m,alpha)
%
abs(cap-5.597545663702324)
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
%
hold on;box on
plot(real(et0),imag(et0),'-k','LineWidth',1.5);
plot(real(et1),imag(et1),'-b','LineWidth',1.5);
%
plot(real(p0(1:2:end)),imag(p0(1:2:end)),'sr','Markerfacecolor','r');
plot(real(p0(2:2:end)),imag(p0(2:2:end)),'dr','Markerfacecolor','r');
%
plot(real(p1(1:2:end)),imag(p1(1:2:end)),'sr','Markerfacecolor','r');
plot(real(p1(2:2:end)),imag(p1(2:2:end)),'dr','Markerfacecolor','r');
%
grid on; grid('minor')
set(gca, 'XMinorTick','on'); set(gca, 'YMinorTick','on')
ax=gca; ax.GridAlpha=0.5; ax.MinorGridAlpha=0.5;
set(gca,'FontSize',22)
set(gca,'LooseInset',get(gca,'TightInset'))
axis equal
axis([-2.2  2.2 -2.2  2.2])
xticks([-2:1:2])
print  -depsc sec-3ex-cirp-m1
%%