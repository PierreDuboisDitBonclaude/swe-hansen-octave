function [direc_log,P, neg]=directivity_view(A,r_zero,lin,ax)

  if nargin < 4 || isempty(ax) || ~ishandle(ax)
    figure;
    ax = axes;
  end

    theta=[A.grid.theta_grid;A.grid.theta_grid(1,:)];
    phi=[A.grid.phi_grid;A.grid.phi_grid(1,:)+2*pi];
    P=squeeze([A.Poynt(:,:,1);A.Poynt(1,:,1)]);
    direc=P/A.P*4*pi*A.r^2;
    direc(direc==0)=realmin;
    direcnorm=direc/max(abs( direc(:) ));
    direc_log=10*log10(direc);

    if min(direc(:))<0
      lin=1;
      neg=1;
    else
      neg=0;
    endif

    if r_zero==0
      a=1;
    elseif r_zero==1
      a=0.5;
    endif

    [X,Y,Z]=sph2cart(phi,pi/2-theta,ones(size(phi))*r_zero+a*direcnorm);

    if lin==0
      surf(ax,X,Y,Z,direc_log);
      if r_zero==0
        zlabel(ax,'Directivity shape r = D / max(|D|)','FontWeight','bold','FontSize',11);
      elseif r_zero==1
        zlabel(ax,'Directivity shape','FontWeight','bold','FontSize',11);
      endif
      cb=colorbar(ax);
      xlabel(cb,'Directivity [dBi]','FontWeight','bold','FontSize',11);
      colormap(ax,viridis);
      caxis(ax,[ 0 max(direc_log(:)) ]);
      ticks=linspace( 0, max(direc_log(:)),6);
      set(cb,'ytick',ticks);
    elseif lin==1
      surf(ax,X,Y,Z,direc);
      if r_zero==0
        zlabel(ax,'Directivity shape r = D / max(|D|)','FontWeight','bold','FontSize',11);
      elseif r_zero==1
        zlabel(ax,'Directivity shape','FontWeight','bold','FontSize',11);
      endif

      cb=colorbar(ax);
      xlabel(cb,'Directivity linear','FontWeight','bold','FontSize',11);
      colormap(ax,viridis);
      caxis(ax,[ min(direc(:)) max(direc(:)) ]);
      #caxis(ax,[ min(a*direc(:)) max(a*direc(:)) ]);
      ticks=linspace( min(direc(:)), max(direc(:)),6);
      set(cb,'ytick',ticks);
    endif
endfunction


