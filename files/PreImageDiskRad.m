function map = PreImageDiskRad (Lc,Lk,ratio,n,tol,Maxiter)
%
%
%
%%
t             =   (0:2*pi/n:2*pi-2*pi/n).'; 
m             =   length(Lc)-1;
thet(1:n,1)   =   pi/2;
for k=2:m+1
    thet(1+(k-1)*n:k*n,1)  =  0;
end
angk     =  angle(Lc);
cent     =  Lc;
radx     = (1-0.5*ratio).*Lk;
rady     =  ratio.*radx;
%%

%%
alpha   =  0;
et(1:n,1)   =   exp(i.*t);
etp(1:n,1)  =   i.*exp(i.*t);
err = inf;
itr = 0;
while (err>tol)
    itr  =itr+1;  
    %%
    for k=2:m+1
        et(1+(k-1)*n:k*n,1)    =  cent(k)+0.5.*exp(i*angk(k)).*(+radx(k).*cos(t)-i*rady(k).*sin(t));
        etp(1+(k-1)*n:k*n,1)   =          0.5.*exp(i*angk(k)).*(-radx(k).*sin(t)-i*rady(k).*cos(t));    
    end
    %%
    A       =  exp(i.*(pi/2-thet)).*(et-alpha);
    for k=1:m+1
        gam((k-1)*n+1:k*n,1)=imag(exp(-i*thet((k-1)*n+1:k*n,1)).*clog(et((k-1)*n+1:k*n,1)-alpha)); 
    end
    %%
    [mun , h ]    =  fbie(et,etp,A,gam,n,5,[],1e-14,200);
    h0            =  sum(h(1:n,1))/n;
    c             =  exp(-h0);
    R             =  h0.*sin(thet)-h;
    %%
    fnet          = (gam+h+i.*mun)./A;
    wn            =  c.*(et-alpha).*exp((et-alpha).*fnet);
    rotwn         =  exp(-i.*R).*wn;
    for k=2:m+1
        Rk(k,1)     =  sum(R((k-1)*n+1:k*n,1))/n; 
        wnL         =  rotwn((k-1)*n+1:k*n,1);
        centk(k,1)  =  exp(i.*Rk(k)).*((max(real(wnL))+min(real(wnL)))/2+i.*(max(imag(wnL))+min(imag(wnL)))/2);     
        radk(k,1)   =  max(real(wnL))-min(real(wnL)); 
    end
    centk;
    radk;
    cent   =  cent-1.0.*(centk-Lc);
    radx   =  radx-(1-0.5*ratio).*(radk-Lk) ;
    angk   =   angle(cent);
    rady   =  ratio.*radx;
    err    = (norm(centk-Lc,1)+norm(radk-Lk,1))/m;
    [itr err]
    error (itr,1) = err;
    itrk  (itr,1) = itr;
    %%
    if itr>=Maxiter
        'No convergence after Maximunm number of iterations'
        break;
    end
end
%%
for k=2:m+1
    et(1+(k-1)*n:k*n,1)    =  cent(k)+0.5.*exp(i*angk(k)).*(+radx(k).*cos(t)-i*rady(k).*sin(t));
    etp(1+(k-1)*n:k*n,1)   =          0.5.*exp(i*angk(k)).*(-radx(k).*sin(t)-i*rady(k).*cos(t));    
end
%%
A       =  exp(i.*(pi/2-thet)).*(et-alpha);
for k=1:m+1
    gam((k-1)*n+1:k*n,1)=imag(exp(-i*thet((k-1)*n+1:k*n,1)).*clog(et((k-1)*n+1:k*n,1)-alpha)); 
end
[mun , h ]    =  fbie(et,etp,A,gam,n,5,[],1e-14,200);
h0            =  sum(h(1:n,1))/n;
c             =  exp(-h0);
fnet          = (gam+h+i.*mun)./A;
zet           =  c.*et.*exp(et.*fnet);
%%
map.zet   =  zet;
map.et    =  et;
map.etp   =  etp;
map.cent  =  cent;
map.radx  =  radx;
map.fnet  =  fnet;
map.c     =  c;
%
end