#Solve_EQS
function Q=solve_eqs(Int_smnc,C,R,phys,grid)
  k=phys.k;
  N=grid.N;
  M=grid.M;
  #Build R_tilde for EQS
  icol=1;
  irow=1;
for col=C
  for row=C
    eqs(icol,irow,:,:)=R_c_gamma(N,col,row,k,R);
    irow+=1;
  endfor
  irow=1;
  icol+=1;
endfor
  eqs_sn=permute(eqs,[4,3,1,2]);
  eqs_sn=permute(repmat(reshape(eqs_sn,[1 size(eqs_sn)]),[M,1,1]),[2 1 3]);

if length(C)==1
  Q   = Int_smnc ./ eqs_sn;
end

##if length(C)==2
##  # Solve EQS for Q_smnc's
##  a = eqs_sn(:,:,1,1);  b = eqs_sn(:,:,1,2);
##  c = eqs_sn(:,:,2,1);  d = eqs_sn(:,:,2,2);
##
##  detA = a.*d - b.*c;
##
##  id1 = reshape(( d - b) ./ detA, size(a,1), 1, size(a,2), 1);
##  id2 = reshape(( a - c) ./ detA, size(a,1), 1, size(a,2), 1);
##
##  Q1 = id1 .* Int_smnc;
##  Q2 = id2 .* Int_smnc;
##
##  Q  = cat(5, Q1, Q2);   % optional: S×M×N×C×2
##endif

