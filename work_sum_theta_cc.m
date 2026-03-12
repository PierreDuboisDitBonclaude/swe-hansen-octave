function out = work_sum_theta_cc(to, grid,phys,fct)

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
  sc=fct.sc;
  sgnM=fct.sgnM;
  z=fct.z;
  R=fct.R;

  [nGrid, mGrid] = meshgrid(0:N, 0:N);
  mask = (mGrid <= nGrid);
  Pcat   = zeros(M, N, blocksize);
  PDxcat = zeros(M, N, blocksize);

  P_mat   = zeros(MN0, MN0);
  P_vec=zeros(1,(N+1)*(N+2)/2);
  P   = zeros(M, N, blocksize);
  [n,m]=meshgrid(1:N,0:N);
  [m_minus, m_plus_0]=NORM(N);
  m_plus=m_plus_0(2:end,:);
  n=n(2:end,:);
  m=m(2:end,:);

  t1  = min(to + blocksize - 1, ts);
  idx2 = to:t1;

  lin = sub2ind([MN0, MN0], mGrid(mask)+1, nGrid(mask)+1);
  l   = nGrid(mask).*(nGrid(mask)+1)/2 + mGrid(mask) + 1;
  idx = [MN0:-1:2, 1:MN0];
  P1  = zeros(MN0, MN0);
  PDx = zeros(MN0, MN0);

            for kk=idx2
                P_vec = gsl_sf_legendre_array(2, MN0, mu(kk), 1);
                P_mat(lin)   = P_vec(l);
                P(:,:,kk-idx2(1)+1)   = P_mat(idx,2:end);
            endfor

            P3_cat=cat( 1, P(MN0:end,:,:) , zeros(1,N,blocksize) );
            dP2 = 1/2 *( ( n - m + 1 ) .* ( n + m ) .* m_minus .*  P3_cat(1:end-2,:,1:length(idx2)) - m_plus .*  P3_cat(3:end,:,1:length(idx2)) );
            dP2a=cat( 1, -m_plus_0(1,:) .* P3_cat(2,:,1:length(idx2)),dP2 );
            P3_=1/2 * reshape( sc(idx2),1,1,length(idx2) ) .* ( ( n - m + 1 ) .* ( n + m ) .* m_minus .*  P3_cat(1:end-2,:,1:length(idx2)) + m_plus .* P3_cat(3:end,:,1:length(idx2)) ) + m.*  reshape(s(idx2),1,1,length(idx2)) .* P3_cat(2:end-1,:,1:length(idx2));
            P3_=cat(1,zeros(1,N,length(idx2)),P3_);
            P3  =  sign(sgnM) .*P3_(idx,:,:);
            f1_a=1i.*P3;
            f2_a=dP2a(idx,:,:);

            f_r=sum(P(:,:,1:length(idx2)).*reshape(abr1,M,N,1).*reshape(squeeze(Q(2,:,:,c)),[M,N,1]),2);
            f_t=sum(f1_a.*reshape(Q1a1,[M,N,1])+f2_a.*reshape(Q2a2,[M,N,1]),2);
            f_p=sum(-f2_a.*reshape(Q1a1,[M,N,1])+f1_a.*reshape(Q2a2,[M,N,1]),2);


  out.sum_r(c,length(R),1:length(idx2),:)=( z.'*squeeze(f_r) ).';
  out.sum_t(c,length(R),1:length(idx2),:)=( z.'*squeeze(f_t) ).';
  out.sum_p(c,length(R),1:length(idx2),:)=( z.'*squeeze(f_p) ).';
  out.idx2 = idx2;
end
