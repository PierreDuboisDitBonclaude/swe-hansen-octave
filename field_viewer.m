function h = field_viewer(orig, recon)


  if ~isstruct(orig) || ~isstruct(recon)
    error('field_viewer3 erwartet zwei Structs: orig und recon');
  end

  % =========================
  % Figure
  % =========================
  h.fig = figure( ...
      'Name','FieldViewer', ...
      'NumberTitle','off', ...
      'Color',[1 1 1], ...
      'Units','normalized', ...
      'Position',[0.04 0.08 0.92 0.84]);

  bg = get(h.fig,'Color');

  btnBg       = [0.84 0.84 0.84];
  btnBg2      = [0.80 0.80 0.80];
  btnBgActive = [0.72 0.72 0.72];

  fsInfo = 12;
  fsForm = 13;

  % =========================
  % Layout
  % =========================
  m = 0.02;
  g = 0.01;

  leftW  = 0.62;
  rightW = 1 - 2*m - g - leftW;

  leftX  = m;
  rightX = leftX + leftW + g;

  h.pLeft  = uipanel('Parent',h.fig,'Units','normalized', ...
                     'Position',[leftX m leftW 1-2*m], ...
                     'BorderType','line','Title','');

  h.pRight = uipanel('Parent',h.fig,'Units','normalized', ...
                     'Position',[rightX m rightW 1-2*m], ...
                     'BorderType','none','Title','');

  % =========================
  % LEFT: Main Axes + Colorbar
  % =========================
  padL = 0.05; padB = 0.08; padT = 0.06; padR = 0.03;
  cbW   = 0.06;
  cbGap = 0.01;

  axW = 1 - padL - padR - cbW - cbGap;
  axH = 1 - padT - padB;

  h.axMain = axes('Parent',h.pLeft,'Units','normalized', ...
                  'Position',[padL padB axW axH]);
  set(h.axMain,'Box','on','LineWidth',1,'FontWeight','bold','FontSize',10);

  imagesc(h.axMain, rand(40,60));
  axis(h.axMain,'xy');
  title(h.axMain,'Main View','FontWeight','bold','FontSize',11);

  h.cbMain = colorbar(h.axMain);
  set(h.cbMain,'FontSize',12,'FontWeight','bold');
  ylabel(h.cbMain,'Directivity [dBi]');

  axpos = get(h.axMain,'Position');
  set(h.cbMain,'Position',[axpos(1)+axpos(3)+cbGap, axpos(2), cbW, axpos(4)]);

  % =========================
  % RIGHT: Controls
  % =========================
  ctlW = 0.52;
  h.pControls = uipanel('Parent',h.pRight,'Units','normalized', ...
                        'Position',[0 0 ctlW 1], ...
                        'BorderType','none','Title','');

  y = 0.995;

  % ============================================================
  % POYNTING:
  % ============================================================
  pTopH = 0.18;

  h.pPTop = uipanel('Parent',h.pControls,'Units','normalized', ...
                    'Position',[0.00 y-pTopH 1.00 pTopH], ...
                    'BorderType','line','Title','');

  h.btnPower = uicontrol('Parent',h.pPTop,'Style','pushbutton', ...
                         'Units','normalized','Position',[0.08 0.62 0.84 0.24], ...
                         'String','3D Directivity View', ...
                         'BackgroundColor',btnBg2, ...
                         'FontWeight','bold','FontSize',9, ...
                         'Callback',@onPoynting);

  h.chkR = uicontrol('Parent',h.pPTop,'Style','checkbox', ...
                     'Units','normalized','Position',[0.08 0.32 0.84 0.18], ...
                     'String','Sphere: r=1 ± 0.5/max(|D|)∙D', ...
                     'BackgroundColor',bg, ...
                     'FontWeight','bold','FontSize',9, ...
                     'Value',0, ...
                     'Callback',@onRChanged);

  h.chkLin = uicontrol('Parent',h.pPTop,'Style','checkbox', ...
                       'Units','normalized','Position',[0.08 0.08 0.84 0.18], ...
                       'String','Color: Log scale', ...
                       'BackgroundColor',bg, ...
                       'FontWeight','bold','FontSize',9, ...
                       'Value',0, ...
                       'Callback',@onLinChanged);

  y = y - pTopH - 0.015;

  % ============================================================
  % POYNTING:
  %
  % ============================================================
  pTextH = 0.23;

  h.pPText = uipanel('Parent',h.pControls,'Units','normalized', ...
                     'Position',[0.00 y-pTextH 1.00 pTextH], ...
                     'BorderType','line','Title','');

  h.axPText = axes('Parent',h.pPText,'Units','normalized', ...
                   'Position',[0.06 0.06 0.92 0.90], ...
                   'Visible','off','HitTest','off');
  set(h.axPText,'XLim',[0 1],'YLim',[0 1]);

  % Info-Zeilen
  h.tF   = text(h.axPText,0,0.95,'','Units','normalized', ...
                'Interpreter','tex','FontSize',fsInfo,'FontWeight','bold', ...
                'VerticalAlignment','top','HorizontalAlignment','left');
  h.tR   = text(h.axPText,0,0.83,'','Units','normalized', ...
                'Interpreter','tex','FontSize',fsInfo,'FontWeight','bold', ...
                'VerticalAlignment','top','HorizontalAlignment','left');
  h.tMax = text(h.axPText,0,0.71,'','Units','normalized', ...
                'Interpreter','tex','FontSize',fsInfo,'FontWeight','bold', ...
                'VerticalAlignment','top','HorizontalAlignment','left');
  h.tP   = text(h.axPText,0,0.59,'','Units','normalized', ...
                'Interpreter','tex','FontSize',fsInfo,'FontWeight','bold', ...
                'VerticalAlignment','top','HorizontalAlignment','left');

  % Formeln/Hinweis (alles in diesem EINEN Block)
  h.tF1  = text(h.axPText,0,0.40,'','Units','normalized', ...
                'Interpreter','tex','FontSize',fsForm, ...
                'VerticalAlignment','top','HorizontalAlignment','left','FontWeight','bold');
  h.tF2  = text(h.axPText,0,0.26,'','Units','normalized', ...
                'Interpreter','tex','FontSize',fsForm, ...
                'VerticalAlignment','top','HorizontalAlignment','left','FontWeight','bold');
  h.tF3  = text(h.axPText,0,0.12,'','Units','normalized', ...
                'Interpreter','tex','FontSize',fsForm, ...
                'VerticalAlignment','top','HorizontalAlignment','left','FontWeight','bold');

  y = y - pTextH - 0.03;

  % =========================
  % Cut Panel
  % =========================
  cutH = 0.275;
  h.pCut = uipanel('Parent',h.pControls,'Units','normalized', ...
                   'Position',[0.00 y-cutH 1.00 cutH], ...
                   'BorderType','line','Title','');

  h.btnCut = uicontrol('Parent',h.pCut,'Style','pushbutton', ...
                       'Units','normalized','Position',[0.08 0.80 0.84 0.16], ...
                       'String','Directivity Cut', ...
                       'BackgroundColor',btnBg2, ...
                       'FontWeight','bold','FontSize',9, ...
                       'Callback',@onCut);

  % Labels
  uicontrol('Parent',h.pCut,'Style','text','Units','normalized', ...
            'Position',[0.10 0.62 0.30 0.10], 'String','Phi', ...
            'BackgroundColor',bg,'FontWeight','bold','FontSize',9, ...
            'HorizontalAlignment','left');

  uicontrol('Parent',h.pCut,'Style','text','Units','normalized', ...
            'Position',[0.56 0.62 0.34 0.10], 'String','Theta-Range', ...
            'BackgroundColor',bg,'FontWeight','bold','FontSize',9, ...
            'HorizontalAlignment','left');

  % Editfelder
  h.edPhi = uicontrol('Parent',h.pCut,'Style','edit','Units','normalized', ...
                      'Position',[0.08 0.52 0.38 0.12], 'String','0', 'FontSize',9, ...
                      'Callback',@onPhiEnter);

  h.edTheta = uicontrol('Parent',h.pCut,'Style','edit','Units','normalized', ...
                        'Position',[0.54 0.52 0.38 0.12], 'String','180', 'FontSize',9, ...
                        'Callback',@onThetaEnter);


  h.txtPhiAct = uicontrol('Parent',h.pCut,'Style','text','Units','normalized', ...
                          'Position',[0.08 0.40 0.38 0.09], ...
                          'String','Phi: -', ...
                          'BackgroundColor',bg,'FontWeight','bold','FontSize',9, ...
                          'HorizontalAlignment','left');

  h.txtThetaAct = uicontrol('Parent',h.pCut,'Style','text','Units','normalized', ...
                            'Position',[0.54 0.40 0.38 0.09], ...
                            'String','Theta: -', ...
                            'BackgroundColor',bg,'FontWeight','bold','FontSize',9, ...
                            'HorizontalAlignment','left');


  h.btnAdd = uicontrol('Parent',h.pCut,'Style','pushbutton','Units','normalized', ...
                       'Position',[0.08 0.10 0.38 0.18], 'String','Add', ...
                       'BackgroundColor',btnBg,'FontWeight','bold','FontSize',9, ...
                       'Callback',@onAddCut);

  h.btnClear = uicontrol('Parent',h.pCut,'Style','pushbutton','Units','normalized', ...
                         'Position',[0.54 0.10 0.38 0.18], 'String','Clear', ...
                         'BackgroundColor',btnBg,'FontWeight','bold','FontSize',9, ...
                         'Callback',@onClearCut);

  y = y - cutH - 0.03;

  % =========================
  % ORIG / RECON
  % =========================
  bH = 0.055;

  h.btnOrig = uicontrol('Parent',h.pControls,'Style','pushbutton', ...
                        'Units','normalized','Position',[0.08 y-bH 0.40 bH], ...
                        'String','ORIG', ...
                        'BackgroundColor',btnBgActive, ...
                        'FontWeight','bold','FontSize',9, ...
                        'Callback',@onOrig);

  h.btnRecon = uicontrol('Parent',h.pControls,'Style','pushbutton', ...
                         'Units','normalized','Position',[0.52 y-bH 0.40 bH], ...
                         'String','RECON', ...
                         'BackgroundColor',btnBg, ...
                         'FontWeight','bold','FontSize',9, ...
                         'Callback',@onRecon);

  y = y - bH - 0.03;

  % =========================
  % Error View
  % =========================
##  h.btnError = uicontrol('Parent',h.pControls,'Style','pushbutton', ...
##                         'Units','normalized','Position',[0.08 y-bH 0.84 bH], ...
##                         'String','Error View', ...
##                         'BackgroundColor',btnBg, ...
##                         'FontWeight','bold','FontSize',9, ...
##                         'Callback',@onError);

  % =========================
  % State
  % =========================
  h.data.orig  = orig;
  h.data.recon = recon;

  h.state.view   = 'poynting';
  h.state.which  = 'orig';
  h.state.r_zero = 0;
  h.state.lin    = 0;

  % CUT: mask-only
  h.state.cutMask = [];
  h.state.cut_phi_in = 0;
  h.state.cut_theta_max = 180;

  guidata(h.fig,h);

  updateMain(); % initial draw

  % ============================================================
  % Callbacks
  % ============================================================
  function onPoynting(~,~)
    hh = guidata(h.fig);
    hh.state.view = 'poynting';
    guidata(h.fig,hh);
    updateMain();
  end

  function onCut(~,~)
    hh = guidata(h.fig);

    cla(hh.axMain,'reset');
    colorbar('off');
    lg = findobj(hh.fig,'Type','legend'); if ~isempty(lg), delete(lg); end

    if strcmp(hh.state.which,'orig'), A = hh.data.orig; else, A = hh.data.recon; end

    if isempty(hh.state.cutMask)
      nPhi = size(A.grid.phi_grid,1);
      hh.state.cutMask = false(nPhi,1);
    end

    hh.state.cut_phi_in = A.grid.phi_grid(1,1);
    hh.state.cut_theta_max = 180;
    set(hh.edPhi,'String',num2str(hh.state.cut_phi_in));
    set(hh.edTheta,'String',num2str(hh.state.cut_theta_max));

    hh.state.cutMask(:) = false;
    hh.state.cutMask(1) = true;

    set(hh.txtPhiAct,'String',sprintf('Phi: %.3g deg', A.grid.phi_grid(1,1)));
    set(hh.txtThetaAct,'String',sprintf('Theta: %.3g deg', hh.state.cut_theta_max));

    hh.state.view = 'cut';
    guidata(h.fig,hh);
    updateMain();
  end

  function onOrig(~,~)
    hh = guidata(h.fig);
    hh.state.which = 'orig';
    set(hh.btnOrig,'BackgroundColor',btnBgActive);
    set(hh.btnRecon,'BackgroundColor',btnBg);
    guidata(h.fig,hh);
    updateMain();
  end

  function onRecon(~,~)
    hh = guidata(h.fig);
    hh.state.which = 'recon';
    set(hh.btnOrig,'BackgroundColor',btnBg);
    set(hh.btnRecon,'BackgroundColor',btnBgActive);
    guidata(h.fig,hh);
    updateMain();
  end

  function onRChanged(src,~)
    hh = guidata(h.fig);
    hh.state.r_zero = double(logical(get(src,'Value')));
    guidata(h.fig,hh);
    updateMain();
  end

  function onLinChanged(src,~)
    hh = guidata(h.fig);
    hh.state.lin = 1-double(logical(get(src,'Value')));
    guidata(h.fig,hh);
    updateMain();
  end

  function onPhiEnter(src,~)
    hh = guidata(h.fig);
    if ~strcmp(hh.state.view,'cut'), return; end

    v = str2double(get(src,'String'));
    if ~isfinite(v), return; end
    hh.state.cut_phi_in = v;

    if strcmp(hh.state.which,'orig'), A = hh.data.orig; else, A = hh.data.recon; end

    phi_vec = A.grid.phi_grid(:,1);
    [~,row] = min(abs(phi_vec/pi*180 - v));
    phi_used = phi_vec(row)/pi*180;

    set(hh.txtPhiAct,'String',sprintf('Phi: %.3g deg', phi_used));

    guidata(h.fig,hh);
  end

  function onThetaEnter(src,~)
    hh = guidata(h.fig);
    if ~strcmp(hh.state.view,'cut'), return; end

    v = str2double(get(src,'String'));
    if ~isfinite(v), return; end
    hh.state.cut_theta_max = v;

    set(hh.txtThetaAct,'String',sprintf('Theta: %.3g deg', hh.state.cut_theta_max));

    guidata(h.fig,hh);
    updateMain();
  end

  function onAddCut(~,~)
    hh = guidata(h.fig);
    if ~strcmp(hh.state.view,'cut'), return; end

    if strcmp(hh.state.which,'orig'), A = hh.data.orig; else, A = hh.data.recon; end

    if isempty(hh.state.cutMask)
      nPhi = size(A.grid.phi_grid,1);
      hh.state.cutMask = false(nPhi,1);
    end

    phi_vec = A.grid.phi_grid(:,1);
    [~,row] = min(abs(phi_vec/pi*180 - hh.state.cut_phi_in));
    phi_used = phi_vec(row)/pi*180;

    hh.state.cutMask(row) = true;

    set(hh.txtPhiAct,'String',sprintf('Phi: %.3g deg', phi_used));
    set(hh.txtThetaAct,'String',sprintf('Theta: %.3g deg', hh.state.cut_theta_max));

    hh.state.view = 'cut';
    guidata(h.fig,hh);
    updateMain();
  end

  function onClearCut(~,~)
    hh = guidata(h.fig);
    if ~strcmp(hh.state.view,'cut'), return; end

    if ~isempty(hh.state.cutMask)
      hh.state.cutMask(:) = false;
    end
    set(hh.txtPhiAct,'String','Phi: -');
    set(hh.txtThetaAct,'String','Theta: -');

    cla(hh.axMain);
    set(hh.tMax,'String','Peak Directivity = -');
    set(hh.tR,'String','r = -');
    set(hh.tP,'String','P_{rad} = -');
    guidata(h.fig,hh);
  end

  % ============================================================
  % Render
  % ============================================================
  function updateMain()
    hh = guidata(h.fig);

    if strcmp(hh.state.which,'orig')
      A = hh.data.orig;
    else
      A = hh.data.recon;
    end

    cla(hh.axMain);

    if strcmp(hh.state.view,'poynting')
      D = [];
      try
        try
          [D, ~, neg] = directivity_view(A, hh.state.r_zero, hh.state.lin, hh.axMain);
        catch
          D = directivity_view(A, hh.state.r_zero, hh.state.lin, hh.axMain);
        end
      catch
      end
      if neg==1
        hh.state.lin = 1;
        set(hh.chkLin,'Value',0);      % Log colouring AUS -> also linear
        set(hh.chkLin,'Visible','off');
      else
        set(hh.chkLin,'Visible','on');
        set(hh.chkLin,'Value', double(~logical(hh.state.lin)));
      endif

      guidata(h.fig,hh);
      title(hh.axMain, sprintf('3D Directivity (%s)', upper(hh.state.which)), ...
            'FontWeight','bold','FontSize',11);

        tighten_axes_to_surface(hh.axMain);



      if isfield(A,'phys') && isfield(A.phys,'f') && isfinite(A.phys.f)
        set(hh.tF,'String',sprintf('f = %.3f GHz', A.phys.f/10^9));
      else
        set(hh.tF,'String','f = -');
      end


      if isfield(A,'r') && isfinite(A.r)
        set(hh.tR,'String',sprintf('r = %.3f m',A.r));
      else
        set(hh.tR,'String','r = -');
      end

      Dok = D(isfinite(D(:)) & imag(D(:))==0 & real(D(:))>0);

      if ~isempty(Dok)
        mx = max(real(Dok));
        set(hh.tMax,'String',sprintf('Peak Directivity = %.2f dBi',mx));
      else
        set(hh.tMax,'String','Peak Directivity = -');
      end

      if isfield(A,'P') && isfinite(A.P)
        set(hh.tP,'String',sprintf('P_{rad} = %.6g W',A.P));
      else
        set(hh.tP,'String','P_{rad} = -');
      end


      f1 = 'D(r,\theta,\phi) = (4\pi r^{2} S_{r}(r,\theta,\phi)) / (\int_{\Omega} S_{r}(r,\theta,\phi) r^{2} d\Omega)';
      f2 = 'P_{rad}(r) = \int_{\Omega} S_{r}(r,\theta,\phi) r^{2} d\Omega';
      f3 = '';
      set(hh.tF1,'String',f1);
      set(hh.tF2,'String',f2);
      set(hh.tF3,'String',f3);

    elseif strcmp(hh.state.view,'cut')
      if isempty(hh.state.cutMask)
        axis(hh.axMain,'off');
        text(hh.axMain,0.02,0.98,'Directivity Cut (empty)', ...
             'Units','normalized','VerticalAlignment','top', ...
             'FontName','monospaced','FontSize',10);
        return;
      end

      phi_vec   = A.grid.phi_grid(:,1);
      theta_vec = A.grid.theta_grid(1,:);
      theta_idx = theta_vec/pi*180 <= hh.state.cut_theta_max;
      if ~any(theta_idx), theta_idx(:) = true; end

      cla(hh.axMain);
      hold(hh.axMain,'on')

      rows = find(hh.state.cutMask);
      peakD = -inf;
      peakTheta = NaN;

      for k = 1:length(rows)
        r = rows(k);
        data = squeeze(A.Poynt(r,:,1));
        direc = data/A.P*4*pi*A.r^2;

        if min(direc(:))<0
          data = direc(:).';
          clin=1;
        else
          direc_log=10*log10(direc);
          data = direc_log(:).';
          clin=0;
        endif
        theta_plot = theta_vec(theta_idx)/pi*180;
        data_plot  = data(theta_idx);

        [localPeak, idxPeak] = max(data_plot);
        localTheta = theta_plot(idxPeak);

        if localPeak > peakD
          peakD = localPeak;
          peakTheta = localTheta;
        end

        plot(hh.axMain, theta_plot, data_plot, 'LineWidth',1.5, ...
             'DisplayName', sprintf('Phi = %.3g deg', phi_vec(r)/pi*180));
        grid on;
      end
      hold(hh.axMain,'off')

      xlabel(hh.axMain,'Theta [deg]','FontWeight','bold','FontSize',11)
      if min(direc(:))<0
        ylabel(hh.axMain,'[1]','FontWeight','bold','FontSize',11)
      else
        ylabel(hh.axMain,'dBi','FontWeight','bold','FontSize',11)
      endif

      lg=legend(hh.axMain,'show');
      set(lg,'FontSize',11,'FontWeight','bold');
      title(hh.axMain, sprintf('Directivity Cut (%s)', upper(hh.state.which)), ...
            'FontWeight','bold','FontSize',12);


      if isfield(A,'phys') && isfield(A.phys,'f') && isfinite(A.phys.f)
        set(hh.tF,'String',sprintf('f = %.3f GHz', A.phys.f/10^9));
      else
        set(hh.tF,'String','f = -');
      end

      if isfield(A,'r') && isfinite(A.r)
        set(hh.tR,'String',sprintf('r = %.3f m',A.r));
      else
        set(hh.tR,'String','r = -');
      end

      if clin==0
        set(hh.tMax,'String',sprintf('Peak Directivity = %.2f dBi (@ %.2f deg)', peakD, peakTheta));
      elseif clin==1
        set(hh.tMax,'String',sprintf('Peak Directivity = %.2f [1] (@ %.2f deg)', peakD, peakTheta));
      else
        set(hh.tMax,'String','Peak Directivity = -');
      end

      if isfield(A,'P') && isfinite(A.P)
        set(hh.tP,'String',sprintf('P_{rad} = %.6g W',A.P));
      else
        set(hh.tP,'String','P_{rad} = -');
      end

      f1 = 'D(r,\theta,\phi) = (4\pi r^{2} S_{r}(r,\theta,\phi)) / (\int_{\Omega} S_{r}(r,\theta,\phi) r^{2} d\Omega)';
      f2 = 'P_{rad}(r) = \int_{\Omega} S_{r}(r,\theta,\phi) r^{2} d\Omega';
      f3 = '';
      set(hh.tF1,'String',f1);
      set(hh.tF2,'String',f2);
      set(hh.tF3,'String',f3);

    else % 'error'
      imagesc(hh.axMain, rand(40,60));
      axis(hh.axMain,'xy');
      title(hh.axMain, sprintf('Error View (%s) [placeholder]', upper(hh.state.which)), ...
            'FontWeight','bold','FontSize',11);
    end

    drawnow;
  end

end


function tighten_axes_to_surface(ax)
  if nargin < 1 || isempty(ax), ax = gca; end

  try
    out=evalc('set(ax,"DataAspectRatio",[1 1 1]);');
  catch err
  end

  s = findobj(ax,'Type','surface');


  if isempty(s), return; end
  s = s(1);

  try
    x = get(s,'XData');
    y = get(s,'YData');
    z = get(s,'ZData');

    xmin = min(x(:)); xmax = max(x(:));
    ymin = min(y(:)); ymax = max(y(:));
    zmin = min(z(:)); zmax = max(z(:));

    if all(isfinite([xmin xmax ymin ymax zmin zmax]))
      padX = 0.02*max(1, xmax-xmin);
      padY = 0.02*max(1, ymax-ymin);
      padZ = 0.02*max(1, zmax-zmin);

      set(ax,'XLim',[xmin-padX, xmax+padX]);
      set(ax,'YLim',[ymin-padY, ymax+padY]);
      set(ax,'ZLim',[zmin-padZ, zmax+padZ]);
    end
  catch
  end
end
