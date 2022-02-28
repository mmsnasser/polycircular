% makefig12.m
% 27-10-2021
clear
addpath ../bie; addpath ../fmm; addpath ../files; addpath ../pcm
%%
rho = @(x,y)(2*asinh(abs(x-y)/(sqrt(1-abs(x)^2)*sqrt(1-abs(y)^2))));
%%
%
n    =  2^11;
t    = [0:2*pi/n:2*pi-2*pi/n].';
et0  =  exp(i.*t); 
et0p =  i.*exp(i.*t); 
%
zv1  = [0.01+0.7i; -0.7+0.01i; -0.01-0.7i;  0.7-0.01i];
zv2  = [0.7+0.01i; -0.01+0.7i; -0.7-0.01i;  0.01-0.7i];
zv3  = (0.7/sqrt(2))*[1+i      ; -1+i      ; -1-i      ;  1-i      ];
zv4  = (0.3/sqrt(2))*[1+i      ; -1+i      ; -1-i      ;  1-i      ];
%
et   =  et0;
etp  =  et0p;
alphav = [];
alpha  = 0;
m      = length(zv1);
%
for k=1:m
    z1   =   zv1(k);
    z2   =   zv2(k);
    z3   =   zv3(k);
    z4   =   zv4(k);
    %
    ver  = [z1  z3  z2  z4];
    dir  = [-1   -1];
    [et1,et1p] = cirarcp3pt(ver,dir,n/2);
    %
    cenz   =  (z1+z2)/2;
    alphav = [alphav;cenz];
    %
    et  = [et  ; et1 ];
    etp = [etp ; et1p];
end
%%
figure(1);
clf
hold on
box on
plot(real(et0),imag(et0),'-k','LineWidth',1.2);
for k=1:m
    J=k*n+1:(k+1)*n;
    crv = et(J);
    plot(real(crv),imag(crv),'-b','LineWidth',1.2);
end
% plot(real(zv1),imag(zv1),'db','Markerfacecolor','b');
% plot(real(alphav),imag(alphav),'dr');
% plot(real(alpha),imag(alpha),'pb');
grid on; grid('minor')
set(gca, 'XMinorTick','on'); set(gca, 'YMinorTick','on')
ax=gca; ax.GridAlpha=0.5; ax.MinorGridAlpha=0.5;
set(gca,'FontSize',22)
axis equal
set(gca,'LooseInset',get(gca,'TightInset'))
axis([-1.1  1.1 -1.1  1.1])
drawnow
print -depsc fig-lens-3
%%
% 
deltv  = [0;1;1;1;1];
% 
format long
cap    = capm(et,etp,alphav,deltv,m,alpha)
capD   = 2*pi/log(1/0.7)
format short g
%