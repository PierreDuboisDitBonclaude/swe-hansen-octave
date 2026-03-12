# ZETA MxL with L=length(Phi)
function z=ZETA(N,phi)
  zeros(2*N+1,length(phi));
  m=(-N:1:N)';
  z=exp(i*m*phi);
end
