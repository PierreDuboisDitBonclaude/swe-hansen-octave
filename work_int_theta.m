function out = work_int_theta(to, fct, phys, grid)
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

       for kk=idx2
          [P_vec, dP_dmu_vec] = gsl_sf_legendre_deriv_array(2, MN0, mu(kk), 1);
          P_mat(lin)   = P_vec(l);
          dP_dmu_mat(lin) = dP_dmu_vec(l);
          P(:,:,kk-idx2(1)+1)   = P_mat(idx,2:end);
          dP_dmu(:,:,kk-idx2(1)+1) = dP_dmu_mat(idx,2:end);
       endfor


      P2    = sum( P(:,:,1:length(idx2)).*reshape(Cr(:,idx2),M,1,length(idx2)) ,3);
      dP2   = sum( dP_dmu(:,:,1:length(idx2)) .* -reshape(s(idx2),1,1,length(idx2)).* reshape(Cp(:,idx2),M,1,length(idx2)) ,3) ;
      P3    = sum( sgnM.* ( P(:,:,1:length(idx2))./reshape(s(idx2),1,1,length(idx2)) ).*reshape(Ct(:,idx2),M,1,length(idx2)) ,3) ;
      dP2b  = sum( dP_dmu(:,:,1:length(idx2)) .* -reshape(s(idx2),1,1,length(idx2)).*reshape(Ct(:,idx2),M,1,length(idx2)) ,3);
      P3b   = sum( sgnM.* ( P(:,:,1:length(idx2))./reshape(s(idx2),1,1,length(idx2)) ).*reshape(Cp(:,idx2),M,1,length(idx2)) ,3);

      out.int_ef1(:,:,1) =ar1.*(1i.*P3 - dP2);
      out.int_ef2(:,:,1) =abr1.*P2 + ar2.* (dP2b + 1i.* P3b);
      out.idx2 = idx2;
end
