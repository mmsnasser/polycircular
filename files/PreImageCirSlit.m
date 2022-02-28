function map = PreImageCirSlitb2(radk,alpk,betk,rr,n,tol,Maxiter)
%
%
%
%%
t        = (0:2*pi/n:2*pi-2*pi/n).'; 
m        =  length(radk);
%%
radk = radk(:); alpk = alpk(:); betk = betk(:);
for k=1:m
    if betk(k)<alpk(k)
        error();
    end
end
%
rk       = (1-0.5*(betk-alpk)/pi).*radk.*(betk-alpk)/pi;
thek     = (alpk+betk)/2;
cent     = (1-0.25*(betk-alpk)/pi).*radk.*exp(i*thek);
rrs      =  rr*(radk/max(radk)).*(2-(betk-alpk)/pi)
%%
for k=1:m
    if abs(cent(k))+rk(k)>0.95
        error('===Choose samll value of radius: ratio===')
    end
end    
%%
et(1:n,1)   =    exp(i*t);
etp(1:n,1)  =  i*exp(i*t);
%
figure(1);
hold on;box on; axis equal
plot(real(et(1:n)),imag(et(1:n)),'k','LineWidth',1.5);
plot(0,0,'pm','LineWidth',1.5);
for k=1:m
    arc  = radk(k).*exp(i.*linspace(alpk(k),betk(k),100));
    plot(real(arc),imag(arc),'g','LineWidth',3);
end
%
err = inf;
itr = 0;
while (err>tol)
    itr  =itr+1;  
    %
    for k=1:m
        Jk = 1+k*n:(k+1)*n;
        et (Jk,1)   =  cent(k)+rk(k)*exp(-i*t);
        etp(Jk,1)   =       -i*rk(k)*exp(-i*t);
    end
    %
    %
    A       =   et;
    gam     =  -log(abs(et));
    %
    [mun,h] = fbie(et,etp,A,gam,n,5,[],1e-14,100);
    %
    fnet    = (gam+h+i.*mun)./A;
    c       =  exp(-mean(h(1:n)));
    zet     =  c*et.*exp(et.*fnet);
    %
figure(1);
for k=1:m
    Jk = 1+k*n:(k+1)*n;
    plot(real(zet(Jk)),imag(zet(Jk)),'-.b','LineWidth',1.5);
end    
for k=1:m
    Jk = 1+k*n:(k+1)*n;
    plot(real(et(Jk)),imag(et(Jk)),'r','LineWidth',1.5);
end    
drawnow    
    %    
    for k=1:m
        Jk = 1+k*n:(k+1)*n;
        Rk(k,1)  =  exp(mean(h(Jk))-mean(h(1:n))); 
        angk     =  carg(zet(Jk));
        ak(k,1)  =  min(angk);
        bk(k,1)  =  max(angk);
        ck(k,1)  = (ak(k)+bk(k))/2;
    end
    %
    Dcent = -(Rk.*exp(i*ck)-radk.*exp(i*thek));
    cent  =  cent + rrs.*Dcent;   
    %
    Drk   = -(Rk.*(bk-ak)-radk.*(betk-alpk))/pi;
    rk    =  rk   + rrs.*Drk;    
      %
%     Dcent = -(Rk-radk)%-((bk-ak)-(betk-alpk))/pi;
%     cent  =  (abs(cent) + rrs.*Dcent).*exp(i*thek);   
    %
%     Drk   =  -0.5*(1-0.5*(betk-alpk)/pi).*(sin(ak-alpk)+sin(betk-bk))
%     rk    =  rk   + rrs.*Drk;    
    %    
    err   = (norm(Rk-radk,1)+norm(alpk-ak,1)+norm(betk-bk,1))/m;
    %
    %
    [itr err]
    errv (itr,1) = err;
    itrk  (itr,1) = itr;
    %
    if itr>=Maxiter
        'No convergence after Maximunm number of iterations'
        break;
    end
end
%%
for k=1:m
    Jk = 1+k*n:(k+1)*n;
    et (Jk,1)   =  cent(k)+rk(k)*exp(-i*t);
    etp(Jk,1)   =       -i*rk(k)*exp(-i*t);
end
%
A       =   et;
gam     =  -log(abs(et));
%
[mun,h] = fbie(et,etp,A,gam,n,5,[],1e-14,100);
%
fnet    = (gam+h+i.*mun)./A;
c       =  exp(-mean(h(1:n)));
zet     =  c*et.*exp(et.*fnet);
%
%%
map.zet   =  zet;
map.et    =  et;
map.etp   =  etp;
map.cent  =  cent;
map.rad   =  rk;
map.fnet  =  fnet;
%%
end