% makefig7c.m
clear
addpath ../bie; addpath ../fmm; addpath ../files
% 
p0   = [5-5i   9+3i  5+7i   9i  -5+7i  -9+5i  -8  -7-5i  -3-8i  -6i  5-5i  ];
d0   = [ 1     1     1      1   -1];
% 
p1   = [-2   -4+i  -6   -4+2i -2];
d1   = [ 1   -1];
% 
p2   = [ 6    4+i   2    2+2i  6];
d2   = [ 1   -1];
% 
p3   = [ 7+5i    5+4i    2+5i    5+6i   7+5i];
d3   = [-1      -1];
% 
p4   = [ 6i   -1+5i   -2+2.5i   -3+4i   -5+5i  -3+7i  6i];
d4   = [ 1     1   -1];
% 
p5   = [ 0   2-i   5-i   -4i  -5-2i  -3-2i 0];
d5   = [ 1    -1   1];
% 
% 
m      =  5;
alphav = [-4+1.5i ; 4+2i ; 4.5+5i ; -2.5+5.5i ;-2i];
alpha  =  2i; 
% 
%%
excap =  21.25386107367032;
nv    = 15.*[2,6,1e1,2e1,5e1,1e2,2e2,5e2,1e3,2e3,4e3,7e3];
res   = [];
for kk=1:length(nv)
    n =   nv(kk)
    %
    [et0,et0p] = cirarcp3pt(p0,d0,n/5);
    [et1,et1p] = cirarcp3pt(p1,d1,n/2);
    [et2,et2p] = cirarcp3pt(p2,d2,n/2);
    [et3,et3p] = cirarcp3pt(p3,d3,n/2);
    [et4,et4p] = cirarcp3pt(p4,d4,n/3);
    [et5,et5p] = cirarcp3pt(p5,d5,n/3);
    %
    et  = [et0  ; et1  ; et2  ; et3  ; et4  ; et5 ];  
    etp = [et0p ; et1p ; et2p ; et3p ; et4p ; et5p];
    %
    cap = capac(et,etp,alphav,m,alpha)
    err = abs(excap-cap)
    res = [res; err];
    %
end
%%
format short g
[nv.' res]
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
%
loglog(nv,res,'-b','LineWidth',2);
hold on;box on
loglog(nv,430*(nv).^(-2.63),'-.k','LineWidth',1.5);
% text(2.3e2,1e-8,'$n^{-2.60}$','FontSize',20)
%
legend('Error','$O(n^{-2.63})$','Location','northeast')
xlabel('$n$')
%
axis square
grid on; grid('minor')
set(gca, 'XMinorTick','on'); set(gca, 'YMinorTick','on')
ax=gca; ax.GridAlpha=0.5; ax.MinorGridAlpha=0.5;
set(gca,'FontSize',20)
set(gca,'LooseInset',get(gca,'TightInset'))
axis([3e1  1e5 1e-12  1e0])
xticks([1e1 1e2 1e3 1e4 1e5])
yticks([1e-12 1e-9 1e-6 1e-3  1e0])
print  -depsc sec-3ex-cirp-m5-err
%%