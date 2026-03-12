function result=get_power2(A,grid)
  if nargin==1
    P=sum(sum(sum(A.*conj(A),1),2),3)*1/2;
    result=P;
    fprintf('\nPower: %d Watt  \n',P);
  elseif nargin>=2 && isstruct(A)
    E=A.E;
    H=A.H;
    Poynt_const=real(cross(E,conj(H),3))/2;
    #Poynt_dyn_comp=cross(E,H,3)/2;
    A.P=sum(Poynt_const(:,:,1)*grid.dphi,1)*grid.wmu*A.r^2;
    A.Poynt=Poynt_const;
    result=A;
    fprintf('\nPower: %d Watt  \n',A.P);
  else
    error('wrong argument')
  endif
endfunction

