
function Int_smnc = Integration(N, phys, grid, field,C,R,name,cores)
  nproc=cores;
  #Integration over phi via FFT
  C1=fftshift(2*pi*ifft(field.E_r,[],1),1);
  C2=fftshift(2*pi*ifft(field.E_theta,[],1),1);
  C3=fftshift(2*pi*ifft(field.E_phi,[],1),1);


  theta=grid.theta;
  mu=grid.mu;
  wmu=grid.wmu;
  ts=grid.ts;

  M=grid.M;
  MN0=grid.MN0;
  N=grid.N;


  k=phys.k;
  eta_0_admit=phys.eta_0_admit;

  Cr=C1.*(wmu.');
  Ct=C2.*(wmu.');
  Cp=C3.*(wmu.');

  ar1=ALPHA(N).*reshape(R1(R,N,C,k),1,N,length(C));
  abr1=ALPHA(N).*reshape(BETA(R,N,k),1,N,1).*reshape(R1(R,N,C,k),1,N,length(C));
  ar2=ALPHA(N).*reshape(R2(R,N,C,k),1,N,length(C));


  sgnM=(-N:N)';
  s=sin(theta);
  sc=cos(theta);


  P2   = zeros(M, N);
  dP2 = zeros(M, N);
  dP2a = zeros(M, N);
  dP2b = zeros(M, N);
  dP2c=zeros(M, N);
  P3 = zeros(M, N);
  P3_ = zeros(M, N);
  P3b = zeros(M, N);
  S1 = zeros(M, N);
  S2 = zeros(M, N);

  [nGrid, mGrid] = meshgrid(0:N, 0:N);
  mask = (mGrid <= nGrid);
  lin = sub2ind([MN0, MN0], mGrid(mask)+1, nGrid(mask)+1);
  l   = nGrid(mask).*(nGrid(mask)+1)/2 + mGrid(mask) + 1;
  idx = [MN0:-1:2, 1:MN0];

  if ( strcmp(name,'block') && strcmp(grid.name,'gl') )
          ## Blocking
          blocksize = 64;

          P_mat   = zeros(MN0, MN0);
          dP_dmu_mat = zeros(MN0, MN0);
          P_vec=zeros(1,(N+1)*(N+2)/2);
          dP_dmu_vec=zeros(1,(N+1)*(N+2)/2);
          P   = zeros(M, N, blocksize);
          dP_dmu = zeros(M, N, blocksize);

          for to = 1:blocksize:ts
              t1  = min(to + blocksize - 1, ts);
              idx2 = to:t1;                  % Theta-Block-Indizes
              for kk=idx2
                  [P_vec, dP_dmu_vec] = gsl_sf_legendre_deriv_array(2, MN0, mu(kk), 1);
                  P_mat(lin)   = P_vec(l);
                  dP_dmu_mat(lin) = dP_dmu_vec(l);
                  P(:,:,kk-idx2(1)+1)   = P_mat(idx,2:end);
                  dP_dmu(:,:,kk-idx2(1)+1) = dP_dmu_mat(idx,2:end);
              endfor

              P2    = sum( P(:,:,1:length(idx2)).*reshape(Cr(:,idx2),M,1,length(idx2)) ,3);  #MxMN0 .*Mx1(1=theta)
              dP2   = sum( dP_dmu(:,:,1:length(idx2)) .* -reshape(s(idx2),1,1,length(idx2)).* reshape(Cp(:,idx2),M,1,length(idx2)) ,3) ;
              P3    = sum( sgnM.* ( P(:,:,1:length(idx2))./reshape(s(idx2),1,1,length(idx2)) ).*reshape(Ct(:,idx2),M,1,length(idx2)) ,3) ;
              dP2b  = sum( dP_dmu(:,:,1:length(idx2)) .* -reshape(s(idx2),1,1,length(idx2)).*reshape(Ct(:,idx2),M,1,length(idx2)) ,3);
              P3b   = sum( sgnM.* ( P(:,:,1:length(idx2))./reshape(s(idx2),1,1,length(idx2)) ).*reshape(Cp(:,idx2),M,1,length(idx2)) ,3);
              S1(:,:) +=ar1.*(1i.*P3 - dP2);
              S2(:,:) +=abr1.*P2 + ar2.*(dP2b + 1i.* P3b);
          endfor

  elseif ( strcmp(name,'para')  && strcmp(grid.name,'gl') )

                pkg load parallel
                warning("off","all");
                more off;
                page_screen_output(0);

                fct.ar1=ar1;
                fct.abr1=abr1;
                fct.ar2=ar2;
                fct.Cr=Cr;
                fct.Ct=Ct;
                fct.Cp=Cp;
                fct.s=s;
                fct.sgnM=sgnM;
                fct.C=C;
                fct.R=R;

                #nproc     = 7;
                fct.blocksize = 16;
                tos = 1:fct.blocksize:ts;
                jobs = num2cell(tos);
                res = parcellfun(nproc, @(to) work_int_theta(to, fct, phys, grid), jobs, "UniformOutput", false);

                for t = 1:numel(res)
                  S1(:,:) += squeeze(res{t}.int_ef1);
                  S2(:,:) += squeeze(res{t}.int_ef2);
                end

  elseif ( strcmp(name,'block') && strcmp(grid.name,'cc') )
          ## Blocking
          blocksize = 64;     % z.B. 8, 16, 32 testen

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

          for to = 1:blocksize:ts
              t1  = min(to + blocksize - 1, ts);
              idx2 = to:t1;
              for kk=idx2
                  P_vec = gsl_sf_legendre_array(2, MN0, mu(kk), 1);
                  P_mat(lin)   = P_vec(l);
                  P(:,:,kk-idx2(1)+1)   = P_mat(idx,2:end);
              endfor
              P3_cat=cat( 1, P(MN0:end,:,:) , zeros(1,N,blocksize) );
              P2    = sum( P(:,:,1:length(idx2)).*reshape(Cr(:,idx2),M,1,length(idx2)) ,3);
              dP2   = 1/2 *( ( n - m + 1 ) .* ( n + m ) .* m_minus .* P3_cat(1:end-2,:,1:length(idx2)) - m_plus .* P3_cat(3:end,:,1:length(idx2)) );
              dP2a  = cat( 1,-m_plus_0(1,:) .* P3_cat(2,:,1:length(idx2)),dP2 );

              dP2b = sum( dP2a(idx,:,:) .* reshape(Cp(:,idx2),M,1,length(idx2)) ,3);
              dP2c = sum( dP2a(idx,:,:) .* reshape(Ct(:,idx2),M,1,length(idx2)) ,3);

              P3_ = 1/2 * reshape( sc(idx2),1,1,length(idx2) ) .* ( ( n - m + 1 ) .* ( n + m ) .* m_minus .* P3_cat(1:end-2,:,1:length(idx2)) + m_plus .* P3_cat(3:end,:,1:length(idx2)) ) + m.* reshape(s(idx2),1,1,length(idx2)) .* P3_cat(2:end-1,:,1:length(idx2));
              P3_ = cat(1,zeros(1,N,length(idx2)),P3_);
              P3  = sum( sign(sgnM) .* P3_(idx,:,:) .* reshape( Ct(:,idx2),M,1,length(idx2) ) ,3);
              P3b = sum( sign(sgnM) .* P3_(idx,:,:) .* reshape( Cp(:,idx2),M,1,length(idx2) ) ,3);

              S1(:,:) +=ar1.*(1i.*P3 - dP2b);
              S2(:,:) +=abr1.*P2 + ar2.*(dP2c + 1i.* P3b);


          endfor
  elseif ( strcmp(name,'para') && strcmp(grid.name,'cc') )
                pkg load parallel
                warning("off","all");
                more off;
                page_screen_output(0);

                fct.ar1=ar1;
                fct.abr1=abr1;
                fct.ar2=ar2;
                fct.Cr=Cr;
                fct.Ct=Ct;
                fct.Cp=Cp;
                fct.s=s;
                fct.sc=sc;
                fct.sgnM=sgnM;
                fct.C=C;
                fct.R=R;

                #nproc     = 7;
                fct.blocksize = 16;
                tos = 1:fct.blocksize:ts;
                jobs = num2cell(tos);
                res = parcellfun(nproc, @(to) work_int_theta_cc(to, fct, phys, grid), jobs, "UniformOutput", false);

                for t = 1:numel(res)
                  S1(:,:) += squeeze(res{t}.int_ef1);
                  S2(:,:) += squeeze(res{t}.int_ef2);
                end
  else
                error('wrong computation method')
  endif

  Int_smnc_u(1,:,:,:)=S1;
  Int_smnc_u(2,:,:,:)=S2;

for mi=-N:N
    int_faktor(mi+N+1)=(-1)^mi;
endfor
int_faktor2=(eta_0_admit^0.5)/k.*int_faktor;
Int_smnc=flip(Int_smnc_u,2).*reshape(int_faktor2,[1,M,1,1]);

end

