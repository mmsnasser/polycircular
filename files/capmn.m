function cap = capmn(et,etp,alphav,deltav,m,alpha)
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
n=length(et)/(m+1); tht=zeros(size(et)); tht(m*n+1:end)=pi/2;
A=exp(-i.*tht).*(et-alpha);
for k=1:m
    for j=1:m+1
        jv = 1+(j-1)*n:j*n;
        gamk{k}(jv,1)=real(exp(-i.*tht(jv)).*clog((et(jv)-alphav(k))./(et(jv)-alpha)));
    end    
    [~,h{k}]=fbie(et,etp,A,gamk{k},n,5,[],1e-14,200);
    for j=1:m+1
        hjk(j,k)=mean(h{k}(1+(j-1)*n:j*n));
    end
end
% Computing the constants a_k  for k=1,2,...,m
mat=hjk; mat(1:m,m+1)=1; mat(m+1,m+1)=0; 
mat(1:m,m+2)=0; mat(m+1,m+2)=-1; 
rhs(1:m,1)=deltav; rhs(m+1,1)=0; 
mat(m+2,1:m)=1; mat(m+2,m+1:m+2)=0; rhs(m+2,1)=0;
x=mat\rhs; a=x(1:m,1); c=x(m+1);
% size(mat)
% cond_numb = cond(mat)
% Computing the capacity
cap  = (2*pi)*sum(deltav(:).*a(:));
% 
end