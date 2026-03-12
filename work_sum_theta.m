function out = work_sum_theta(to, grid,phys,fct)
  warning("off","all");
  more off;
  page_screen_output(0);

  theta=grid.theta;
  mu=grid.mu;
  ts=grid.ts;
  M=fct.M;
  MN0=fct.MN0;
  N=fct.N;

  k=phys.k;
  eta_0_admit=phys.eta_0_admit;

  blocksize=fct.blocksize;
  Q1a1=fct.Q1a1;
  abr1=fct.abr1;
  Q2a2=fct.Q2a2;
  Q=fct.Q;
  c=fct.c;
  s=fct.s;
  sgnM=fct.sgnM;
  z=fct.z;
  R=fct.R;

  [nGrid, mGrid] = meshgrid(0:N, 0:N);
  mask = (mGrid <= nGrid);
  Pcat   = zeros(M, N, blocksize);
  PDxcat = zeros(M, N, blocksize);

  t1  = min(to + blocksize - 1, ts);
  idx2 = to:t1;

  lin = sub2ind([MN0, MN0], mGrid(mask)+1, nGrid(mask)+1);
  l   = nGrid(mask).*(nGrid(mask)+1)/2 + mGrid(mask) + 1;
  idx = [MN0:-1:2, 1:MN0];
  P1   = zeros(MN0, MN0);
  PDx = zeros(MN0, MN0);

  for kk=idx2
    [Pvec, dPdx_vec] = gsl_sf_legendre_deriv_array(2, N, mu(kk), 1);
    P1(lin)   = Pvec(l);
    PDx(lin) = dPdx_vec(l);
    Pcat(:,:,kk-idx2(1)+1)   = P1(idx,2:end);
    PDxcat(:,:,kk-idx2(1)+1) = PDx(idx,2:end);
  endfor

  f1_a=(1i.*sgnM.*Pcat(:,:,1:length(idx2))./reshape(s(idx2),1,1,length(idx2)));
  f2_a=(PDxcat(:,:,1:length(idx2)).*-reshape(s(idx2),1,1,length(idx2)));
  f_r=sum(Pcat(:,:,1:length(idx2)).*reshape(abr1,M,N,1).*reshape(squeeze(Q(2,:,:,c)),[M,N,1]),2);
  f_t=sum(f1_a.*reshape(Q1a1,[M,N,1])+f2_a.*reshape(Q2a2,[M,N,1]),2);
  f_p=sum(-f2_a.*reshape(Q1a1,[M,N,1])+f1_a.*reshape(Q2a2,[M,N,1]),2);

  out.sum_r(c,length(R),1:length(idx2),:)=( z.'*squeeze(f_r) ).';
  out.sum_t(c,length(R),1:length(idx2),:)=( z.'*squeeze(f_t) ).';
  out.sum_p(c,length(R),1:length(idx2),:)=( z.'*squeeze(f_p) ).';
  out.idx2 = idx2;
end
