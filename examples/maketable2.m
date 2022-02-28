% maketable2
% 02-06-2021
clear;
addpath ../bie; addpath ../fmm; addpath ../files; addpath ../pcm
%% 
av   =  [0.1, 0.5, 0.1+i*0.3, -0.2+i*0.5, -0.3-i*0.5].';
Ta   =  @(z,a)((z-a)./(1-conj(a).*z));
Tap  =  @(z,a)((1-abs(a)^2)./((1-conj(a).*z).^2));
%
n    =  2^12;
t    = [0:2*pi/n:2*pi-2*pi/n].';
et0  =  exp(i.*t); 
et0p =  i.*exp(i.*t); 
vera = [-0.2   0.6i  0.2   0.1i -0.2];
d1   = [-1   +1];
[eta,etap] = cirarcp3pt(vera,d1,n/2);
%
deltv = [0;1];
et   = [et0 ;  eta  ];
etp  = [et0p ; etap ];
alphav0 =  0.3i;
alpha0  = -0.3i;
cap0    = capm(et,etp,alphav0,deltv,1,alpha0);
%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on
box on
plot(real(et0),imag(et0),'-k','LineWidth',2);
plot(real(eta),imag(eta),'-k','LineWidth',2);
plot(real(vera),imag(vera),'sk','Markerfacecolor','k');
grid on; grid('minor')
set(gca, 'XMinorTick','on'); set(gca, 'YMinorTick','on')
ax=gca; ax.GridAlpha=0.5; ax.MinorGridAlpha=0.5;
set(gca,'FontSize',22)
axis equal
set(gca,'LooseInset',get(gca,'TightInset'))
xticks([-1:0.5:1])
yticks([-1:0.5:1])
axis([-1.1  1.1 -1.1  1.1])
drawnow
%%    
for k=1:length(av)
    a    =  av(k);
    et1  =  Ta(eta,a);
    et1p =  Tap(eta,a).*etap;
    %
    ver1  =  Ta(vera,a);
    %
    et   = [et0 ;  et1  ];
    etp  = [et0p ; et1p ];
    alphav =  Ta(alphav0,a);
    alpha  =  Ta(alpha0,a);
    %
    ff=figure;
    set(groot,'defaultAxesTickLabelInterpreter','latex');  
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');
    hold on;box on;
    plot(real(et0),imag(et0),'-k','LineWidth',2);
    plot(real(et1),imag(et1),'-k','LineWidth',2);
    plot(real(ver1),imag(ver1),'sk','Markerfacecolor','k');
%     plot(real(alphav),imag(alphav),'or');
%     plot(real(alpha),imag(alpha),'pr');
    grid on; grid('minor')
    set(gca, 'XMinorTick','on'); set(gca, 'YMinorTick','on')
    ax=gca; ax.GridAlpha=0.5; ax.MinorGridAlpha=0.5;
    set(gca,'FontSize',22)
    axis equal
    set(gca,'LooseInset',get(gca,'TightInset'))
    xticks([-1:0.5:1])
    yticks([-1:0.5:1])
    axis([-1.1  1.1 -1.1  1.1])
    drawnow
    %
    deltv = [0;1];
    cap(k,1) = capm(et,etp,alphav,deltv,1,alpha);
end
%%
format short g
[0;av]
format long g
jj = [1:length(av)]';
res = [  cap0 NaN];
res = [ res ;    cap abs(cap-cap0)];
res
%%