# R1 NxCxL mit L=length(r)

function r1=R1(R,N,C,K)
  r1=zeros(N,length(C),length(R));
  [j,n,h1,h2]=sph_all_bessel(N, K*R);
  A=[j;n;h1;h2].';
  r1(:,1:length(C),length(R))=A(2:end,C);
end

