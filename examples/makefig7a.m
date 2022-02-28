% makefig7a.m
clear
addpath ../bie; addpath ../fmm; addpath ../files
% 
p0   = [ 1+i  2i -1+i -2 -1-i -2i 1-i  2 1+i];
d0   = [ 1    1   1    1];
%
p1   = [ 0.5+0.5i  0.4  0.5-0.5i -0.4i -0.5-0.5i -0.4 -0.5+0.5i  0.4i  0.5+0.5i];
d1   = [ 1         1    1         1];
%
excap =  5.597545663702171; 
%%
nv  = [10,100,500,1000,5000,10000,50000,100000];
res = [];
for kk=1:length(nv)
    n =   nv(kk)
    [et0,et0p] = cirarcp3pt(p0,d0,n/2);
    [et1,et1p] = cirarcp3pt(p1,d1,n/2);
    %
    et  = [et0  ; et1  ];  
    etp = [et0p ; et1p ];
    alphav = [ 0 ];
    alpha = i; 
    m =1;
    %
    cap = capac(et,etp,alphav,m,alpha)
    err = abs(excap-cap)
    res = [res; err];
    %
end
%%
format short g
[nv.' res]
% tab =[     10      0.33222
%           100    0.0010516
%           500   1.2763e-05
%          1000   1.9156e-06
%          5000   2.3866e-08
%         10000   3.6215e-09
%         50000   4.5307e-11
%         1e+05   6.1684e-12];
%%
f = figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
%
loglog(nv,res,'-b','LineWidth',2);
hold on;box on
%
loglog(nv,320*(nv).^(-2.73),'-.k','LineWidth',1.5);
% loglog(nv,(nv).^(-2.73),'-.k','LineWidth',1.5);
% text(1e2,2e-8,'$n^{-2.73}$','FontSize',20)
legend('Error','$O(n^{-2.73})$','Location','northeast')
%
xlabel('$n$')
%
axis square
grid on; grid('minor')
set(gca, 'XMinorTick','on'); set(gca, 'YMinorTick','on')
ax=gca; ax.GridAlpha=0.5; ax.MinorGridAlpha=0.5;
set(gca,'FontSize',20)
set(gca,'LooseInset',get(gca,'TightInset'))
axis([1e1  1e5 1e-12  1e0])
xticks([1e1 1e2 1e3 1e4 1e5])
yticks([1e-12 1e-9 1e-6 1e-3  1e0])

print  -depsc sec-3ex-cirp-m1-err
%%