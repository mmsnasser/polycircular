function f = halfslitmap(Lc,Lk,thetk,r,n,alphah,tloerance,maxiter )
%
% Nasser, June 11, 2019
% Given an unbounded multiply connected domain Omega of connectivity m+1
% consists of m slits in the lower half-plan where
% Lc: is a vector of the center points of these slits
% Lk: is a vector of the length of these slits
% thetk: is a vector of the angles between these slits and the positive
% real line.
% 
% This function computes: zet , zetp , cent , wet; i.e., a preimage 
% bounded multiply connected circular domain G border by m pseudo-ellipses 
% interior to the unit circle:
% zet: the parametruzation of the circles. 
% zetp: the derivative of zet.
% cent: the center of the circles
% wzet: the boundary values of the conformal mapping from G onto Omega,
% 
% where:
% n: the number of discretization points in each boundary component of G
% r: the ration of the length of the minor axis to the major axis
% tolerance: is the tolerance of the iterative method
% maxiter: is the maximum number of iterations for the iterative method
%
%
% This code of computing the preimage domain G and the conformal mapping
% is based on the iterative method presented in Section 4 in the paper:
% M. Nasser and C. Green, A fast numerical method for ideal fluid flow 
% in domains with multiple stirrers, Nonlinearity 31 (2018) 815-837.
% 
%%
t             =   (0:2*pi/n:2*pi-2*pi/n).'; 
%
Psi       =  @(z)(1i.*(2i./(i-z)-1)); %i(i+z)/(i-z)
Psip      =  @(z)((-2./((i-z).^2)));
Psiv      =  @(z)(i+2./(i+z));
Psivp     =  @(z)(-2./((i+z).^2));
%
Lc = Lc(:);  Lk = Lk(:);  thetk = thetk(:); 
Lc = [0;Lc]; Lk = [0;Lk]; thetk = [0;thetk];
%
m             =   length(thetk)-1;
thet(1:n,1)   =   0;
for k=2:m+1
    thet(1+(k-1)*n:k*n,1)  =  thetk(k);
end
cent     =  Lc;
radx     = (1-0.5*r).*Lk;
rady     =  r.*radx;
%
alpha       =   Psiv(alphah)
et(1:n,1)   =   exp(i.*t);
etp(1:n,1)  =   i.*exp(i.*t);
err = inf;
itr = 0;
% 
while (err>tloerance)
    itr  =itr+1;  
    %
    for k=2:m+1
        xt(1+(k-1)*n:k*n,1)    =  cent(k)+0.5.*exp(i*thetk(k)).*(+radx(k).*cos(t)-i*rady(k).*sin(t));
        xtp(1+(k-1)*n:k*n,1)   =          0.5.*exp(i*thetk(k)).*(-radx(k).*sin(t)-i*rady(k).*cos(t));    
    end
    %
    for k=2:m+1
        et(1+(k-1)*n:k*n,1)    =  Psiv(xt(1+(k-1)*n:k*n,1));
        etp(1+(k-1)*n:k*n,1)   =  Psivp(xt(1+(k-1)*n:k*n,1)).*xtp(1+(k-1)*n:k*n,1);
    end
    %    
    A       =  exp(i.*(pi/2-thet)).*(et-alpha);
    gam(1:n,1)   =   0;
    for k=2:m+1
        gam((k-1)*n+1:k*n,1)=imag(exp(-i*thet((k-1)*n+1:k*n,1)).*Psi(et((k-1)*n+1:k*n,1))); 
    end
    %
    [mun , h ]    =  fbie(et,etp,A,gam,n,5,[],1e-14,100);
    h0            =  sum(h(1:n,1))/n;
    %
    fnet          = (gam+h+i.*mun)./A;
    zet           =  Psi(et)+(et-alpha).*fnet+i*h0;
    rotwn         =  exp(-i.*thet).*zet;
    for k=2:m+1
        wnL         =  rotwn((k-1)*n+1:k*n,1);
        centk(k,1)  =  exp(i.*thetk(k)).*((max(real(wnL))+min(real(wnL)))/2+i.*(max(imag(wnL))+min(imag(wnL)))/2);     
        radk(k,1)   =  max(real(wnL))-min(real(wnL)); 
    end
    cent   =  cent-1.0.*(centk-Lc);
    radx   =  radx-(1-0.5*r).*(radk-Lk) ;
    rady   =  r.*radx;
    err    = (norm(centk-Lc,1)+norm(radk-Lk,1))/(m);
    [itr err]
    error (itr,1) = err;
    itrk  (itr,1) = itr;
    %
    if itr>=maxiter
        'No convergence after Maximunm number of iterations'
        break;
    end
end
%
for k=2:m+1
    xt(1+(k-1)*n:k*n,1)    =  cent(k)+0.5.*exp(i*thetk(k)).*(+radx(k).*cos(t)-i*rady(k).*sin(t));
    xtp(1+(k-1)*n:k*n,1)   =          0.5.*exp(i*thetk(k)).*(-radx(k).*sin(t)-i*rady(k).*cos(t));    
end
for k=2:m+1
    et(1+(k-1)*n:k*n,1)    =  Psiv(xt(1+(k-1)*n:k*n,1));
    etp(1+(k-1)*n:k*n,1)   =  Psivp(xt(1+(k-1)*n:k*n,1)).*xtp(1+(k-1)*n:k*n,1);
    centc(k,1)             =  Psiv(cent(k));
end
A       =  exp(i.*(pi/2-thet)).*(et-alpha);
gam(1:n,1)   =   0;
for k=2:m+1
    gam((k-1)*n+1:k*n,1)=imag(exp(-i*thet((k-1)*n+1:k*n,1)).*Psi(et((k-1)*n+1:k*n,1))); 
end
[mun , h ]    =  fbie(et,etp,A,gam,n,5,[],1e-14,100);
h0            =  sum(h(1:n,1))/n;
fnet          = (gam+h+i.*mun)./A;
zet           =  Psi(et)+(et-alpha).*fnet+i*h0;
%%
f.zet    =  zet;
f.cent   =  centc;
f.et     =  et;
f.etp    =  etp;
f.alpha  =  alpha;
end