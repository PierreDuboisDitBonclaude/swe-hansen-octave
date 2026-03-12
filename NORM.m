
function [Mminus1,Mplus1]=NORM(N)
  A=zeros(N+1,N+1);
  NO_mat=zeros(N+1,N+1);
  NO_mat2=zeros(N+2,N+1);
  Mplus1_=zeros(N+1,N+1);
  Mminus1_=zeros(N+1,N+1);

  [Ng,Mg]=meshgrid(0:N,0:N);

  mask=Mg<=Ng;
  n=Ng(mask);
  m=Mg(mask);

  r_minus = 1./sqrt((n-m+1).*(n+m));
  r_plus  = sqrt((n+m+1).*(n-m));
  Mminus1_(mask)=r_minus;
  Mplus1_(mask)=r_plus;
  Mminus1=triu(Mminus1_)(2:end,2:end);
  Mplus1=triu(Mplus1_,1)(:,2:end);
  Mplus1(1,:)=sqrt(Ng(1,2:end).*(Ng(1,2:end)+1));
endfunction

