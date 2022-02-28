function cap = cap2s(et,etp,alphav,alpha)
% Compute the capacity of the condensers (G,F)
% 
% Input:
% 1,2) et, etp: parametrization of the boundary and its first derivative
% 3) alphav=[alphav(1),...,alphav(mp)]: alphav(j) is an auxiliary point
% interior to \Gamma_j
% 4) m: the number of the closed sets E_k
% 8) alpha: for bounded G, alpha is an auxiliary point in G
% 
% Output:
% cap (the capacity of the condensers (G,F)).
%
% Computing the constants \h_{j,k} for j=1,2,...,m
% 
n=length(et)/3; tht=zeros(size(et)); tht(2*n+1:3*n)=pi/2;
A=exp(-i.*tht).*(et-alpha);
for k=1:2
    for j=1:3
        jv = 1+(j-1)*n:j*n;
        gamk{k}(jv,1)=real(exp(-i.*tht(jv)).*clog((et(jv)-alphav(k))./(et(jv)-alpha)));
    end    
    [~,h{k}]=fbie(et,etp,A,gamk{k},n,5,[],1e-14,200);
    for j=1:3
        hv(j,k)=mean(h{k}(1+(j-1)*n:j*n));
    end
end
% Computing the constants a_k  for k=1,2,...,m
a1 = 1/((hv(1,1)+hv(2,2))-(hv(1,2)+hv(2,1)));
% Computing the capacity
cap  = (2*pi)*a1;
% 
end