
# R2 NxCxL mit L=kength(r)

function r2=R2(R,N,C,K)
  lc=length(C);
  lr=length(R);
  r2=zeros(N,lc,lr);

  r1=zeros(N+1,length(C),length(R));
  [j,n,h1,h2]=sph_all_bessel(N+1, K*R);
  A=[j;n;h1;h2].';
  N1=(1:N).';

  r2=( 1./(2*N1+1).*( (N1+1).*A(N1,:,:) -N1.*A(N1+2,:,:) ) )(:,C);
end

