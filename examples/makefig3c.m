% makefig3c.m
clear
addpath ../bie; addpath ../fmm; addpath ../files
% 
n =   3*5*2^9
p0   = [5-5i   9+3i  5+7i   9i  -5+7i  -9+5i  -8  -7-5i  -3-8i  -6i  5-5i  ];
d0   = [ 1     1     1      1   -1];
[et0,et0p] = cirarcp3pt(p0,d0,n/5);
% 
p1   = [-2   -4+i  -6   -4+2i -2];
d1   = [ 1   -1];
[et1,et1p] = cirarcp3pt(p1,d1,n/2);
% 
p2   = [ 6    4+i   2    2+2i  6];
d2   = [ 1   -1];
[et2,et2p] = cirarcp3pt(p2,d2,n/2);
% 
p3   = [ 7+5i    5+4i    2+5i    5+6i   7+5i];
d3   = [-1      -1];
[et3,et3p] = cirarcp3pt(p3,d3,n/2);
% 
p4   = [ 6i   -1+5i   -2+2.5i   -3+4i   -5+5i  -3+7i  6i];
d4   = [ 1     1   -1];
[et4,et4p] = cirarcp3pt(p4,d4,n/3);
% 
p5   = [ 0   2-i   5-i   -4i  -5-2i  -3-2i 0];
d5   = [ 1    -1   1];
[et5,et5p] = cirarcp3pt(p5,d5,n/3);
% 
et  = [et0  ; et1  ; et2  ; et3  ; et4  ; et5 ];  
etp = [et0p ; et1p ; et2p ; et3p ; et4p ; et5p];
% 
m      =  5;
alphav = [-4+1.5i ; 4+2i ; 4.5+5i ; -2.5+5.5i ;-2i];
alpha  =  2i; 
% 
format long 
cap = capac(et,etp,alphav,m,alpha)
abs(cap-21.25386107367032)
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
plot(real(et3),imag(et3),'-b','LineWidth',1.5);
plot(real(et4),imag(et4),'-b','LineWidth',1.5);
plot(real(et5),imag(et5),'-b','LineWidth',1.5);
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
plot(real(p3(1:2:end)),imag(p3(1:2:end)),'sr','Markerfacecolor','r');
plot(real(p3(2:2:end)),imag(p3(2:2:end)),'dr','Markerfacecolor','r');
%
plot(real(p4(1:2:end)),imag(p4(1:2:end)),'sr','Markerfacecolor','r');
plot(real(p4(2:2:end)),imag(p4(2:2:end)),'dr','Markerfacecolor','r');
%
plot(real(p5(1:2:end)),imag(p5(1:2:end)),'sr','Markerfacecolor','r');
plot(real(p5(2:2:end)),imag(p5(2:2:end)),'dr','Markerfacecolor','r');
%
grid on; grid('minor')
set(gca, 'XMinorTick','on'); set(gca, 'YMinorTick','on')
ax=gca; ax.GridAlpha=0.5; ax.MinorGridAlpha=0.5;
set(gca,'FontSize',22)
set(gca,'LooseInset',get(gca,'TightInset'))
axis equal
axis([-10  10  -8.5  9.5])
xticks([-9:3:9])
yticks([-8:2:8])
print  -depsc sec-3ex-cirp-m5
%%