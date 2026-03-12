function field=get_field2(grid,phys,R,C,Q,name,nproc)

  k=phys.k;
  phi_value=grid.phi.';
  N=size(Q)(3);                           #
  M=size(Q)(2);                           #
  ts=grid.ts;
  ps=grid.ps;
  mu=grid.mu;
  theta_value=grid.theta.';
  sgnM=(-N:N).';
  eta_0_admit=phys.eta_0_admit;
  MN0=N+1;                                #


  z=ZETA(N,phi_value);
  a=ALPHA(N);
  b=BETA(R,N,k);
  r1=R1(R,N,C,k);
  r2=R2(R,N,C,k);
  ab=a.*b.';

  N1=N+1;
  [nGrid, mGrid] = meshgrid(0:N, 0:N);
  mask = (mGrid <= nGrid);
  lin = sub2ind([N1, N1], mGrid(mask)+1, nGrid(mask)+1);
  l   = nGrid(mask).*(nGrid(mask)+1)/2 + mGrid(mask) + 1;
  idx = [N1:-1:2, 1:N1];               % 1×M

  P   = zeros(N1, N1);
  PDx = zeros(N1, N1);
  P2   = zeros(M, N1);
  PDx2 = zeros(M, N1);
  P3 = zeros(M, N1);

  sum_r=zeros(length(C),length(R),ts,ps);
  sum_t=zeros(length(C),length(R),ts,ps);
  sum_p=zeros(length(C),length(R),ts,ps);

  x=mu;
  s=sin(theta_value);
  sc=cos(theta_value);
  ls=length(theta_value);

  if ( strcmp(name,'block') && strcmp(grid.name,'cc') )

    for c=1:length(C)
        ar1=a.*reshape(r1(:,c),1,N);
        abr1=ab.*reshape(r1(:,c),1,N);
        ar2=a.*reshape(r2(:,c),1,N);
        Q1a1=squeeze(Q(1,:,:,c)).*reshape(ar1,M,N);
        Q2a2=squeeze(Q(2,:,:,c)).*reshape(ar2,M,N);

        blocksize = 64;     % z.B. 8, 16, 32 testen

        P_mat   = zeros(MN0, MN0);
        P_vec=zeros(1,(N+1)*(N+2)/2);
        P   = zeros(M, N, blocksize);
        [n,m]=meshgrid(1:N,0:N);
        [m_minus, m_plus_0]=NORM(N);
        m_plus=m_plus_0(2:end,:);
        n=n(2:end,:);
        m=m(2:end,:);

        for to = 1:blocksize:ts
            t1  = min(to + blocksize - 1, ts);
            idx2 = to:t1;

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

            sum_r(c,length(R),idx2,:)=( z.'*squeeze(f_r) ).';
            sum_t(c,length(R),idx2,:)=( z.'*squeeze(f_t) ).';
            sum_p(c,length(R),idx2,:)=( z.'*squeeze(f_p) ).';
        endfor

    endfor

  endif

  if ( strcmp(name,'para') && strcmp(grid.name,'cc') )

    for c=1:length(C)
      ar1=a.*reshape(r1(:,c),1,N);
      abr1=ab.*reshape(r1(:,c),1,N);
      ar2=a.*reshape(r2(:,c),1,N);
      Q1a1=squeeze(Q(1,:,:,c)).*reshape(ar1,M,N);
      Q2a2=squeeze(Q(2,:,:,c)).*reshape(ar2,M,N);

      pkg load parallel

      warning("off","all");
      more off;
      page_screen_output(0);

      #nproc     = 7
      blocksize = 16;

      fct.Q1a1=Q1a1;
      fct.abr1=abr1;
      fct.Q2a2=Q2a2;
      fct.Q=Q;
      fct.c=c;
      fct.s=s;
      fct.sc=sc;
      fct.sgnM=sgnM;
      fct.z=z;
      fct.R=R;
      fct.M=M;
      fct.N=N;
      fct.MN0=MN0;

      tos = 1:blocksize:ts;
      fct.blocksize=blocksize;

      jobs = num2cell(tos);
      res = parcellfun(nproc, @(to) work_sum_theta_cc(to,grid,phys,fct), jobs, "UniformOutput", false);
      sum_r = zeros(ts, ps);
      sum_t = zeros(ts, ps);
      sum_p = zeros(ts, ps);

      for k = 1:numel(res)
        idx2 = res{k}.idx2;
        sum_r(idx2, :) = squeeze(res{k}.sum_r);
        sum_t(idx2, :) = squeeze(res{k}.sum_t);
        sum_p(idx2, :) = squeeze(res{k}.sum_p);
      end

    endfor

  endif

  if ( strcmp(name,'block') && strcmp(grid.name,'gl') )

    for c=1:length(C)
      ar1=a.*reshape(r1(:,c),1,N);
      abr1=ab.*reshape(r1(:,c),1,N);
      ar2=a.*reshape(r2(:,c),1,N);
      Q1a1=squeeze(Q(1,:,:,c)).*reshape(ar1,M,N);
      Q2a2=squeeze(Q(2,:,:,c)).*reshape(ar2,M,N);

      blocksize = 64;
      Pcat   = zeros(M, N, blocksize);
      PDxcat = zeros(M, N, blocksize);
            for to = 1:blocksize:ts
              t1  = min(to + blocksize - 1, ts);
              idx2 = to:t1;

              for kk=idx2
                [Pvec, dPdx_vec] = gsl_sf_legendre_deriv_array(2, N, x(kk), 1);
                P(lin)   = Pvec(l);
                PDx(lin) = dPdx_vec(l);
                Pcat(:,:,kk-idx2(1)+1)   = P(idx,2:end);
                PDxcat(:,:,kk-idx2(1)+1) = PDx(idx,2:end);
              endfor

              f1_a=(1i.*sgnM.*Pcat(:,:,1:length(idx2))./reshape(s(idx2),1,1,length(idx2)));
              f2_a=(PDxcat(:,:,1:length(idx2)).*-reshape(s(idx2),1,1,length(idx2)));

              f_r=sum(Pcat(:,:,1:length(idx2)).*reshape(abr1,M,N,1).*reshape(squeeze(Q(2,:,:,c)),[M,N,1]),2);
              f_t=sum(f1_a.*reshape(Q1a1,[M,N,1])+f2_a.*reshape(Q2a2,[M,N,1]),2);
              f_p=sum(-f2_a.*reshape(Q1a1,[M,N,1])+f1_a.*reshape(Q2a2,[M,N,1]),2);

              sum_r(c,length(R),idx2,:)=( z.'*squeeze(f_r) ).';
              sum_t(c,length(R),idx2,:)=( z.'*squeeze(f_t) ).';
              sum_p(c,length(R),idx2,:)=( z.'*squeeze(f_p) ).';
            end
    end
  endif
  if ( strcmp(name,'para') && strcmp(grid.name,'gl') )
    for c=1:length(C)
      ar1=a.*reshape(r1(:,c),1,N);
      abr1=ab.*reshape(r1(:,c),1,N);
      ar2=a.*reshape(r2(:,c),1,N);
      Q1a1=squeeze(Q(1,:,:,c)).*reshape(ar1,M,N);
      Q2a2=squeeze(Q(2,:,:,c)).*reshape(ar2,M,N);

      pkg load parallel

      warning("off","all");
      more off;
      page_screen_output(0);

      #nproc     = 7
      blocksize = 16;

      fct.Q1a1=Q1a1;
      fct.abr1=abr1;
      fct.Q2a2=Q2a2;
      fct.Q=Q;
      fct.c=c;
      fct.s=s;
      fct.sgnM=sgnM;
      fct.z=z;
      fct.R=R;
      fct.M=M;
      fct.N=N;
      fct.MN0=MN0;

      tos = 1:blocksize:ts;
      fct.blocksize=blocksize;

      jobs = num2cell(tos);
      res = parcellfun(nproc, @(to) work_sum_theta(to,grid,phys,fct), jobs, "UniformOutput", false);
      sum_r = zeros(ts, ps);
      sum_t = zeros(ts, ps);
      sum_p = zeros(ts, ps);

      for k = 1:numel(res)
        idx2 = res{k}.idx2;
        sum_r(idx2, :) = squeeze(res{k}.sum_r);
        sum_t(idx2, :) = squeeze(res{k}.sum_t);
        sum_p(idx2, :) = squeeze(res{k}.sum_p);
      end
    endfor
  endif
      field.sum_r=squeeze(sum_r);
      field.sum_t=squeeze(sum_t);
      field.sum_p=squeeze(sum_p);
      field.grid=grid;
      field.phys=phys;
endfunction

