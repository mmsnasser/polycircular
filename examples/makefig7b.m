% makefig7b.m
clear
addpath ../bie; addpath ../fmm; addpath ../files
% 
n =   2^12
p0   = [-8   -4i    8    8i   -8];
d0   = [ 1    1];
p1   = [-2   -4+i  -6   -4+2i -2];
d1   = [ 1   -1];
p2   = [ 6    4+i   2    2+2i  6];
d2   = [ 1   -1];
alphav = [-4+1.5i ; 4+2i ];
alpha = 4i; m =2;
%
%
excap =  10.95831564100293;
nv    = [10,100,500,1000,5000,10000,50000,100000];
res   = [];
for kk=1:length(nv)
    n =   nv(kk)
    %
    [et0,et0p] = cirarcp3pt(p0,d0,n/2);
    [et1,et1p] = cirarcp3pt(p1,d1,n/2);
    [et2,et2p] = cirarcp3pt(p2,d2,n/2);
    %
    et  = [et0  ; et1  ; et2 ];  
    etp = [et0p ; et1p ; et2p];
    %
    cap = capac(et,etp,alphav,m,alpha)
    err = abs(excap-cap)
    res = [res; err];
    %
end
%%
f = figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
%
loglog(nv,res,'-b','LineWidth',2);
hold on;box on
%
loglog(nv,470*(nv).^(-2.65),'-.k','LineWidth',1.5);
legend('Error','$O(n^{-2.65})$','Location','northeast')
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
%
print  -depsc sec-3ex-cirp-m2-err
%%