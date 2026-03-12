function [h bg] = mode_viewer(Q)
% MODE_VIEWER5(Q)


  h.fig = figure( ...
    'Name','ModeViewer', ...
    'NumberTitle','off', ...
    'Color',[1 1 1], ...
    'Units','normalized', ...
    'Position',[0.05 0.08 0.90 0.84]);

  % -------------------------
  % Global layout
  % -------------------------
  m  = 0.03;
  g  = 0.02;
  rightW = 0.12;

  leftW  = 1 - 2*m - g - rightW;
  leftX  = m;
  rightX = leftX + leftW + g;

  leftY = m;
  leftH = 1 - 2*m;

  h.pLeft  = uipanel('Parent',h.fig,'Units','normalized', ...
                     'Position',[leftX leftY leftW leftH], ...
                     'BorderType','none');

  h.pRight = uipanel('Parent',h.fig,'Units','normalized', ...
                     'Position',[rightX leftY rightW leftH], ...
                     'BorderType','line','Title','');

  % -------------------------
  % Left inner layout
  % -------------------------
  pm = 0.02;
  innerX = pm; innerY = pm;
  innerW = 1 - 2*pm;
  innerH = 1 - 2*pm;

  rowGap = 0.03;

  baseTopFrac = (0.58 * 0.70) * 1.05;
  delta  = 0.05;
  topFrac = baseTopFrac + delta;
  botFrac = (1 - rowGap/innerH) - topFrac;

  topH = topFrac * innerH;
  botH = botFrac * innerH;

  topY = innerY + botH + rowGap;
  botY = innerY;

  colGap = 0.04;
  colW   = (innerW - colGap) / 2;

  % -------------------------
  % TOP PANELS
  % -------------------------
  h.pTop1 = uipanel('Parent',h.pLeft,'Units','normalized', ...
                    'Position',[innerX              topY colW topH], ...
                    'BorderType','line','Title','');
  h.pTop2 = uipanel('Parent',h.pLeft,'Units','normalized', ...
                    'Position',[innerX+colW+colGap  topY colW topH], ...
                    'BorderType','line','Title','');

  h.axPlot1 = axes('Parent',h.pTop1,'Units','normalized', ...
                   'Position',[0.05 0.06 0.90 0.84]);
  h.axPlot2 = axes('Parent',h.pTop2,'Units','normalized', ...
                   'Position',[0.05 0.06 0.90 0.84]);

  set(h.axPlot1,'Box','on','LineWidth',1,'FontWeight','bold','FontSize',9);
  set(h.axPlot2,'Box','on','LineWidth',1,'FontWeight','bold','FontSize',9);

  % -------------------------
  % BOTTOM PANELS
  % -------------------------
  h.pBars1 = uipanel('Parent',h.pLeft,'Units','normalized', ...
                     'Position',[innerX              botY colW botH], ...
                     'BorderType','line','Title','');
  h.pBars2 = uipanel('Parent',h.pLeft,'Units','normalized', ...
                     'Position',[innerX+colW+colGap  botY colW botH], ...
                     'BorderType','line','Title','');

  % -------------------------
  % RIGHT SIDE:
  % -------------------------
  rInnerX = pm; rInnerY = pm;
  rInnerW = 1 - 2*pm;
  rInnerH = 1 - 2*pm;

  ctrlH = 0.28;
  tblH  = rInnerH - ctrlH;

  h.pTbl = uipanel('Parent',h.pRight,'Units','normalized', ...
                   'Position',[rInnerX rInnerY rInnerW tblH], ...
                   'BorderType','line','Title','');

  h.pCtrl = uipanel('Parent',h.pRight,'Units','normalized', ...
                    'Position',[rInnerX rInnerY+tblH rInnerW ctrlH], ...
                    'BorderType','none','Title','');


  colnames = { 'Qsmn', 'P[W]', '%', 'Σ P[W]', 'Σ %' };

  h.tbl = uitable('Parent',h.pTbl,'Units','normalized', ...
                  'Position',[0.03 0.03 0.94 0.94], ...
                  'ColumnName',colnames, ...
                  'RowName',[], ...
                  'ColumnEditable',[false false false false false], ...
                  'ColumnWidth',{80 70 55 70 55}, ...
                  'Data',cell(0,5));

  set(h.tbl, 'FontSize',8)

  % -------------------------
  % Threshold Inputs
  % -------------------------
  bg = get(h.fig,'Color');
  x0 = 0.06; w0 = 0.88;


  infoH   = 0.16;
  rowH    = 0.30;
  gapH    = 0.06;
  topPad  = 0.04;

  yTop = 1 - topPad;


  h.pInfoDummy = uipanel('Parent',h.pCtrl,'Units','normalized', ...
                         'Position',[x0 yTop-infoH w0 infoH], ...
                         'BorderType','line','Title','');

  h.axInfoDummy = axes('Parent',h.pInfoDummy, ...
                       'Units','normalized', ...
                       'Position',[0.05 0.15 0.90 0.70], ...
                       'Visible','off');

  set(h.axInfoDummy,'XLim',[0 1],'YLim',[0 1]);
  P=1/2*sum( (Q.*conj(Q))(:) );
  h.txtInfoDummy = text(h.axInfoDummy,0,0.5, ...
                        sprintf('P = 1/2 ∑ |Q_{smn}|^2 = %.5g W',P), ...
                        'Interpreter','tex', ...
                        'Units','normalized', ...
                        'FontSize',12, ...
                        'FontWeight','bold', ...
                        'HorizontalAlignment','left', ...
                        'VerticalAlignment','middle');

  yTop = yTop - infoH - 0.04;

  labels = {'SMN','SMN'};
  texts  = {' Lower Threshold',' Upper Threshold'};
  edits = cell(1,2);

  for k=1:2
      p = uipanel('Parent',h.pCtrl,'Units','normalized', ...
                  'Position',[x0 yTop-rowH w0 rowH], ...
                  'BorderType','none');

      if k == 1
          qy  = 0.52;
          ly  = 0.40;
          ty  = 0.52;
          ey  = 0.00;
      else
          qy  = 0.46;
          ly  = 0.34;
          ty  = 0.46;
          ey  = 0.00;
      end

      uicontrol('Parent',p,'Style','text','Units','normalized',...
          'Position',[0.00 qy 0.12 0.32],'String','Q',...
          'HorizontalAlignment','left','FontWeight','bold','FontSize',9,...
          'BackgroundColor',bg);

      uicontrol('Parent',p,'Style','text','Units','normalized',...
          'Position',[0.10 ly 0.22 0.32],'String',labels{k},...
          'HorizontalAlignment','left','FontWeight','bold','FontSize',7,...
          'BackgroundColor',bg);

      uicontrol('Parent',p,'Style','text','Units','normalized',...
          'Position',[0.30 ty 0.70 0.32],'String',texts{k},...
          'HorizontalAlignment','left','FontWeight','bold','FontSize',8,...
          'BackgroundColor',bg);

      edits{k} = uicontrol('Parent',p,'Style','edit','Units','normalized',...
          'Position',[0.00 ey 1.00 0.32],'FontSize',9);

      yTop = yTop - rowH - gapH;
  end

  h.e_q_lo = edits{1};
  h.e_q_hi = edits{2};


  set(h.e_q_lo,'String', num2str(0));

  P = real(0.5 * Q .* conj(Q));
  P = P(isfinite(P));
  if isempty(P), pmax = 1; else, pmax = max(P(:)); end
  set(h.e_q_hi,'String', num2str(pmax,17));


  set(h.e_q_lo,'Callback',@(src,evt) update_plots_cb2(h, Q));
  set(h.e_q_hi,'Callback',@(src,evt) update_plots_cb2(h, Q));

  % Optional: initial refresh
  update_plots_cb2(h, Q);

end
