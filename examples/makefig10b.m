% makefig10b.m
% 27-10-2021
clear
addpath ../bie; addpath ../fmm; addpath ../files; addpath ../pcm
%%
rho = @(x,y)(2*asinh(abs(x-y)/(sqrt(1-abs(x)^2)*sqrt(1-abs(y)^2))));
%%
nv   = [10,100,500,1000,5000,10000,50000,100000]';
res  = [];
capv = [];
for kk=1:length(nv)
    n =   nv(kk)
    %
    t    = [0:2*pi/n:2*pi-2*pi/n].';
    et0  =  exp(i.*t); 
    et0p =  i.*exp(i.*t); 
    %
    r    =  4/5;
    s    =  3/10;
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
    cap    = capm(et,etp,alphav,deltv,1,alpha);
    %
    capE = 10.15585205509004;
    err  = abs(cap-capE);
    %
    capv = [capv;cap];
    res  = [res;err];
end
%%
format long g
[nv capv]
res
%
f = figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
%
loglog(nv,res,'-b','LineWidth',2);
hold on;box on
%
loglog(nv,1300*(nv).^(-2.95),'-.k','LineWidth',1.5);
% loglog(nv,(nv).^(-2.73),'-.k','LineWidth',1.5);
% text(1e2,2e-8,'$n^{-2.73}$','FontSize',20)
legend('Error','$O(n^{-2.95})$','Location','northeast')
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
print  -depsc sec-lens-1-err
print  -dpdf  sec-lens-1-err
%%