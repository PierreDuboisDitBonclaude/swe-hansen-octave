function [Zmin, Zmax]=mode_chart4(Q,s,c,th_low,th_high,a,ax,chartTitle)

  if (ndims(Q)==4)
    Q = squeeze(Q(:,:,:,c));
  endif

  if strcmp(a,'pow')
    P_Q = real(0.5*Q.*conj(Q));
    q   = squeeze(P_Q(s,:,:));
  elseif strcmp(a,'phase')
    error('not implemented yet')
  endif

cmap = [1 1 1; jet(256); 0.6 0.6 0.6; 0.15 0.15 0.15];
N = rows(cmap);

maxth=max(P_Q(s,P_Q(s,:)<=th_high));
minth=min(P_Q(s,P_Q(s,:)>=th_low));
if maxth<=minth
  maxth=minth+eps;
end


for i=linspace(0.8,1.2,1000)
  zmax=maxth+i*2*(maxth-minth)/256;
  zmin=minth-i*(maxth-minth)/256;
  T = octave_cmap_bins([zmin zmax], N);
  if (minth>T(2,2) && minth<T(2,3)) && (maxth>T(end-2,2) && maxth<T(end-2,3))
    brk=1;
    break;
  endif
endfor

masklow=T(1,2)+(T(1,3)-T(1,2))/2;
masktop=T(end,2)+(T(end,3)-T(end,2))/2;
maskhigh=T(end-1,2)+(T(end-1,3)-T(end-1,2))/2;

##########################

  q2 = q';
  [ny, nx] = size(q2);

  x = (1:nx) - ceil(nx/2);   % M
  y = 1:ny;                  % N

  [X, Y] = meshgrid(x, y);
  valid = (abs(X) <= Y);

  B = q2;
  B(B>th_high)   = maskhigh;
  B(~valid) = masktop;
  B(B<th_low)   = masklow;

  % ----- draw -----
  imagesc(ax, x, y, B);

  colormap(ax,cmap);
  caxis(ax,[zmin, zmax]);

  % ----- styling -----
  set(ax,"XAxisLocation","bottom");
  box(ax,"on");
  grid(ax,"off");
  set(ax,"LineWidth",1.0,"FontWeight","bold","FontSize",9);

  % ----- title -----
  title(ax, chartTitle, "FontWeight","bold", "FontSize", 12, "Interpreter","tex");

  % ----- labels -----
  ylabel(ax, "N", "FontWeight","bold", "FontSize", 10);
  xlabel(ax, "M", "FontWeight","bold", "FontSize", 10);

  % ----- colorbar (thin + closer) -----
  cb = colorbar(ax);

  title(cb,"P[W]","FontWeight","bold","FontSize",10);

  ylim(cb,[max(th_low,0) min( th_high,max(B(B<=th_high)) )]);

  ticks = linspace(max(th_low,0),min( th_high,max(B(B<=th_high)) ),6);
  set(cb,'ytick',ticks);

  labels = cell(size(ticks));
  for k = 1:numel(ticks)
      labels{k} = sprintf('%.2e', ticks(k));
  end
  set(cb,'yticklabel', labels);

  set(cb,"Units","normalized");
  p = get(cb,"Position");
  p(1) = p(1) - 0.03;
  p(3) = p(3) * 0.65;
  set(cb,"Position",p);


  pos = get(ax,"Position");
  pos(2) = pos(2) + 0.03;
  pos(4) = pos(4) - 0.03;
  set(ax,"Position",pos);

  drawnow();
end
