function out = work_int_theta_cc(to, fct, phys, grid)
        warning("off","all");
        more off;
        page_screen_output(0);

  theta=grid.theta;
  mu=grid.mu;
  wmu=grid.wmu;
  ts=grid.ts;

  M=grid.M;
  MN0=grid.MN0;
  N=grid.N;
  k=phys.k;

  Cr=fct.Cr;
  Ct=fct.Ct;
  Cp=fct.Cp;
  ar1=fct.ar1;
  abr1=fct.abr1;
  ar2=fct.ar2;
  sgnM=fct.sgnM;
  s=fct.s;
  sc=fct.sc;
  blocksize=fct.blocksize;


[nGrid, mGrid] = meshgrid(0:N, 0:N);
mask = (mGrid <= nGrid);


  t1  = min(to + blocksize - 1, ts);
  idx2 = to:t1;
  lin = sub2ind([MN0, MN0], mGrid(mask)+1, nGrid(mask)+1);
  l   = nGrid(mask).*(nGrid(mask)+1)/2 + mGrid(mask) + 1;
  idx = [MN0:-1:2, 1:MN0];

          P2   = zeros(M, N);
          dP2 = zeros(M, N);
          dP2b = zeros(M, N);
          P3 = zeros(M, N);
          P3b = zeros(M, N);

          P_mat   = zeros(MN0, MN0);
          dP_dmu_mat = zeros(MN0, MN0);
          P_vec=zeros(1,(N+1)*(N+2)/2);
          dP_dmu_vec=zeros(1,(N+1)*(N+2)/2);
          P   = zeros(M, N, blocksize);
          dP_dmu = zeros(M, N, blocksize);

          [n,m]=meshgrid(1:N,0:N);
          [m_minus, m_plus_0]=NORM(N);
          m_plus=m_plus_0(2:end,:);
          n=n(2:end,:);
          m=m(2:end,:);

          for kk=idx2
              P_vec = gsl_sf_legendre_array(2, MN0, mu(kk), 1);
              P_mat(lin)   = P_vec(l);
              P(:,:,kk-idx2(1)+1)   = P_mat(idx,2:end);
          endfor

          P3_cat=cat( 1, P(MN0:end,:,:) , zeros(1,N,blocksize) );

          P2    = sum( P(:,:,1:length(idx2)).*reshape(Cr(:,idx2),M,1,length(idx2)) ,3);
          dP2 = 1/2 *( ( n - m + 1 ) .* ( n + m ) .* m_minus .* P3_cat(1:end-2,:,1:length(idx2)) - m_plus .* P3_cat(3:end,:,1:length(idx2)) );
          dP2a=cat( 1,-m_plus_0(1,:) .* P3_cat(2,:,1:length(idx2)),dP2 );

          dP2b = sum( dP2a(idx,:,:) .* reshape(Cp(:,idx2),M,1,length(idx2)) ,3);
          dP2c = sum( dP2a(idx,:,:) .* reshape(Ct(:,idx2),M,1,length(idx2)) ,3);

          P3_=1/2 * reshape( sc(idx2),1,1,length(idx2) ) .* ( ( n - m + 1 ) .* ( n + m ) .* m_minus .* P3_cat(1:end-2,:,1:length(idx2)) + m_plus .* P3_cat(3:end,:,1:length(idx2)) ) + m .* reshape(s(idx2),1,1,length(idx2)) .* P3_cat(2:end-1,:,1:length(idx2));
          P3_=cat(1,zeros(1,N,length(idx2)),P3_);
          P3  = sum( sign(sgnM) .* P3_(idx,:,:) .* reshape( Ct(:,idx2),M,1,length(idx2) ) ,3);
          P3b = sum( sign(sgnM) .* P3_(idx,:,:) .* reshape( Cp(:,idx2),M,1,length(idx2) ) ,3);

      out.int_ef1(:,:,1) =ar1.*(1i.*P3 - dP2b);
      out.int_ef2(:,:,1) =abr1.*P2 + ar2.* (dP2c + 1i.* P3b);
      out.idx2 = idx2;
end
