
# BETA NxL L=length(r)
function b=BETA(R,N,K)
  b=zeros(N,length(R));
##  [N1,~]=meshgrid(0:N,0:N);
##  N2=cat(1,N1,N1)(2:end,:)
  N1=1:N;
  b=(N1.*(N1+1)./(K*R)).';
end
