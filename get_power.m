function [P,Poynt_const,Poynt_dyn_comp]=get_power(A,grid,sw)
  #gets power of Q or E-H-field
  if nargin==1
    P=sum(sum(sum(A.*conj(A),1),2),3)*1/2;
  endif

  if nargin>=2
    E=0;
    H=0;
    P_=0;
    P=0;
    E=A.E;
    H=A.H;
    Poynt_const=real(cross(E,conj(H),3))/2;
    Poynt_dyn_comp=cross(E,H,3)/2;    # only part of dynamic Poynting vector formulation...needs time to be complete....
    P=sum(Poynt_const(:,:,1)*grid.dphi,1)*grid.wmu*A.r^2;
  endif

  if nargin==2
    fprintf('\nPower: %d Watt  \n',P);
  endif

endfunction

