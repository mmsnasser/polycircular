% makefig10a.m
% 27-10-2021
clear
addpath ../bie; addpath ../fmm; addpath ../files; addpath ../pcm
%%
rho = @(x,y)(2*asinh(abs(x-y)/(sqrt(1-abs(x)^2)*sqrt(1-abs(y)^2))));
%%
%
n    =  2^12;
t    = [0:2*pi/n:2*pi-2*pi/n].';
et0  =  exp(i.*t); 
et0p =  i.*exp(i.*t); 
%
r    =  0.8;
s    =  0.3;
%
ver  = [-r   s*i  r   -s*i  ];
dir  = [-1   -1];
[et1,et1p] = cirarcp3pt(ver,dir,n/2);
%
et  = [et0  ; et1 ];
etp = [et0p ; et1p];
%
deltv  = [0;1];
alphav =  0.0i;
alpha  =  (1+s)*i/2;
% 
cap  = capm(et,etp,alphav,deltv,1,alpha);
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
clf
hold on
box on
plot(real(et0),imag(et0),'-k','LineWidth',1.2);
plot(real(et1),imag(et1),'-b','LineWidth',1.2);
plot(real(ver),imag(ver),'db','Markerfacecolor','b');
% plot(real(alphav),imag(alphav),'dr');
% plot(real(alpha),imag(alpha),'pb');
grid on; grid('minor')
set(gca, 'XMinorTick','on'); set(gca, 'YMinorTick','on')
ax=gca; ax.GridAlpha=0.5; ax.MinorGridAlpha=0.5;
set(gca,'FontSize',22)
axis equal
xticks([-1:0.5:1])
yticks([-1:0.5:1])
axis([-1.1 1.1 -1.1 1.1])
set(gca,'FontSize',20)
set(gca,'LooseInset',get(gca,'TightInset'))
print -depsc fig-lens-1
drawnow
%%
[c,R] =  my3Pts(r,s*i,-r);
d     = -imag(c);
%
if r>s
    syms rr
    f         = @(rr)((16*rr./((1-rr.^2).^2)).*asin((R^2-rr.^2-d^2)./(2*d*rr)));
    iit       =  double(vpaintegral(f,rr,s,r));
elseif r==s
    iit  = 0;
end
hA    =  4*pi*(s.^2./(1-s.^2))+iit;
%
% compute the       hyperbolic perimeter       of the Lens   
P  = (2*pi/n)*sum(2.*abs(et1p)./(1-abs(et1).^2));
LB   =  4*pi/(exp(2*mu(tanh(P/4)))-1);
%
%%
