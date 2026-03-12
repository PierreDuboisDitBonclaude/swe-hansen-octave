
function R_c_gamma_s_n=R_c_gamma(N,c,gamma,K,R)
  n=(1:N)';
  R_c_gamma_s_n=zeros(N,2);
  R_c_gamma_s_n(:,1)=squeeze( R1(R,N,c,K).*R1(R,N,gamma,K) );
  R_c_gamma_s_n(:,2)=squeeze( R2(R,N,c,K).*R2(R,N,gamma,K)+(n.*(n+1).*R1(R,N,c,K).*R1(R,N,gamma,K)./(K*R)^2) );
endfunction

