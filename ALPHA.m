# ALPHA MxN
function [a,N2, M2]=ALPHA(N)
  a=zeros(2*N+1,N);
  [N1,M1]=meshgrid(1:N,0:N);
  N2=cat(1,N1,N1)(2:end,:);
  M2=cat(1,-flip(M1(2:end,:),1),M1);

  a=1/((2*pi)^0.5)*1./((N2.*(N2+1)).^0.5).*(-(M2)./abs(M2)).^M2;
  a(N+1,:)=(1/((2*pi)^0.5)*1./((N2.*(N2+1)).^0.5))(N+1,:);

endfunction

