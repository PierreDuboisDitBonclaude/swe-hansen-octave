function field=make_field(Q,R,C,phys,grids,name,cores)

  E=get_field2(grids,phys,R,C,Q,name,cores);
  H=get_field2(grids,phys,R,C,flip(Q,1),name,cores);

  field.E_r=( E.sum_r.*phys.k/((phys.eta_0_admit)^0.5) ).';
  field.E_theta=( E.sum_p.*phys.k/((phys.eta_0_admit)^0.5) ).';
  field.E_phi=( E.sum_p.*phys.k/((phys.eta_0_admit)^0.5) ).';

  field.H_r=H.sum_r.'*phys.k*phys.eta_0_admit^0.5*-i;
  field.H_theta=H.sum_t.'*phys.k*phys.eta_0_admit^0.5*-i;
  field.H_phi=H.sum_p.'*phys.k*phys.eta_0_admit^0.5*-i;

  field.E=cat(3,E.sum_r.',E.sum_t.',E.sum_p.')*phys.k/((phys.eta_0_admit)^0.5);
  field.H=cat(3,H.sum_r.',H.sum_t.',H.sum_p.')*phys.k*phys.eta_0_admit^0.5*-i;

  field.r=R;
  field.grid=E.grid;
  field.phys=phys;
endfunction

