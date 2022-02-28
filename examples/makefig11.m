% makefig11.m
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
sv   = [linspace(0.05,r,31)].';
%%
for kk=1:length(sv)
    s    =  sv(kk);
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
    cap(kk,1)  = capm(et,etp,alphav,deltv,1,alpha);
    %
    figure(1);
    clf
    hold on
    box on
    plot(real(et0),imag(et0),'-k','LineWidth',1.2);
    plot(real(et1),imag(et1),'-k','LineWidth',1.2);
    plot(real(ver),imag(ver),'pk','Markerfacecolor','k');
    plot(real(alphav),imag(alphav),'dr');
    plot(real(alpha),imag(alpha),'pb');
    grid on; grid('minor')
    set(gca, 'XMinorTick','on'); set(gca, 'YMinorTick','on')
    ax=gca; ax.GridAlpha=0.5; ax.MinorGridAlpha=0.5;
    set(gca,'FontSize',22)
    axis equal
    set(gca,'LooseInset',get(gca,'TightInset'))
    axis([-1.1  1.1 -1.1  1.1])
    drawnow
    %
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
    hA(kk,1)  =  4*pi*(s.^2./(1-s.^2))+iit;
    %
    % compute the       hyperbolic perimeter       of the Lens   
    p  = (2*pi/n)*sum(2.*abs(et1p)./(1-abs(et1).^2));
    ps(kk,1) =  p;
    %
end
%%
capD = 2*pi/log(1/r)
capS = 2*pi/mu(2*r/(1+r^2))
figure;
hold on
box on
plot(sv,sv-sv+capD,':k','LineWidth',1.5);
plot(sv,sv-sv+capS,'--k','LineWidth',1.5);
plot(sv,cap,'-b','LineWidth',1.5);
% plot(sv,hA,'-r','LineWidth',1.2);
% plot(sv,LB,'-b','LineWidth',1.2);
% legend({'hyp-area','LB'},'Location','northwest','Interpreter','latex')
xlabel('$s$','Interpreter','latex')
ylabel('Capacity','Interpreter','latex')
grid on; grid('minor')
set(gca, 'XMinorTick','on'); set(gca, 'YMinorTick','on')
ax=gca; ax.GridAlpha=0.5; ax.MinorGridAlpha=0.5;
set(gca,'FontSize',22)
set(gca,'LooseInset',get(gca,'TightInset'))
axis([0  0.8 6  30])
print -depsc fig-lens-2
%%
capD = 2*pi/log(1/r)
capS = 2*pi/mu(2*r/(1+r^2))
figure;
hold on
box on
for k=1:length(ps)
    LB(k,1) = 2*pi/mu(tanh(ps(k)/4));
end
UB = 2*pi./log(sqrt(1+(2*pi./ps).^2)+2*pi./ps);
plot(ps,cap,'-b','LineWidth',1.5);
plot(ps,UB,'-.k','LineWidth',1.5);
plot(ps,LB,'-.r','LineWidth',1.5);
plot(ps,ps-ps+capD,':k','LineWidth',1.5);
plot(ps,ps-ps+capS,'--','color','k','LineWidth',1.5);
legend({'${\rm cap}(D,E)$','UB','LB'},'Location','northwest','Interpreter','latex')
xlabel('$L=$hyp-perim$(E)$','Interpreter','latex')
% ylabel('Capacity','Interpreter','latex')
grid on; grid('minor')
set(gca, 'XMinorTick','on'); set(gca, 'YMinorTick','on')
ax=gca; ax.GridAlpha=0.5; ax.MinorGridAlpha=0.5;
set(gca,'FontSize',22)
set(gca,'LooseInset',get(gca,'TightInset'))
axis([8.75 28 5  30])
print -depsc fig-lens-2hp
%%