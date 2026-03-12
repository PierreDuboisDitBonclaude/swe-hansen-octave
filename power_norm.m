function  Pn_field=power_norm(field,Pn)
  if isstruct(field)
    P_old=get_power(field,field.grid,1);
    Pn_field=field;
    Pnc=sqrt(Pn/P_old);

    Pn_field.E=Pn_field.E*Pnc;
    Pn_field.H=Pn_field.H*Pnc;

    Pn_field.E_r=Pn_field.E_r*Pnc;
    Pn_field.E_theta=Pn_field.E_theta*Pnc;
    Pn_field.E_phi=Pn_field.E_phi*Pnc;

    Pn_field.H_r= Pn_field.H_r*Pnc;
    Pn_field.H_theta= Pn_field.H_theta*Pnc;
    Pn_field.H_phi= Pn_field.H_phi*Pnc;

    [Pn_field.P , Pn_field.Poynt]=get_power(Pn_field,field.grid,1);
    Pn_field.Pnc=Pnc;

    fprintf('\n Normalised Power:       %d Watt\n',Pn_field.P);
    fprintf(  ' Normalisation Constant: %d \n',Pnc);
  else
    P_old=get_power(field);
    Pnc=sqrt(Pn/P_old);
    Pn_field=field*Pnc;
    P=get_power(Pn_field);
    fprintf('\n Normalised Power:       %d Watt\n',P);
    fprintf(  ' Normalisation Constant: %d \n',Pnc);
  endif

endfunction

