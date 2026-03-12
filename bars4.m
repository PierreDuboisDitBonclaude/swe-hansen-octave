function bars4(pnl, Q2, barsPerRow, barTitle)
% bars4(pnl, Q2, barsPerRow, barTitle)


  if nargin < 4, barTitle = ''; end


  ok = isappdata(pnl,'axP') && isappdata(pnl,'axPh') && isappdata(pnl,'axTitle') && ...
       isgraphics(getappdata(pnl,'axP')) && isgraphics(getappdata(pnl,'axPh')) && isgraphics(getappdata(pnl,'axTitle'));

  if ok
    axP     = getappdata(pnl,'axP');
    axPh    = getappdata(pnl,'axPh');
    axTitle = getappdata(pnl,'axTitle');
  else
    delete(get(pnl,'Children'));


    titleH = 0.08;
    padL = 0.06; padR = 0.05;
    padT = 0.02;
    padB = 0.10;
    gap  = 0.10;

    axW = 1 - padL - padR;
    usableH = 1 - titleH - padT - padB;
    axH = (usableH - gap) / 2;


    axP  = axes('Parent',pnl,'Units','normalized', ...
                'Position',[padL padB + axH + gap  axW axH]);
    axPh = axes('Parent',pnl,'Units','normalized', ...
                'Position',[padL padB              axW axH]);

    setappdata(pnl,'axP',axP);
    setappdata(pnl,'axPh',axPh);


    axTitle = axes('Parent',pnl,'Units','normalized', ...
                   'Position',[0 0 1 1], ...
                   'Visible','off', ...
                   'HitTest','off', ...
                   'PickableParts','none');
    setappdata(pnl,'axTitle',axTitle);
  end


  axes(axTitle);
  cla(axTitle);


  text(axTitle, 0.5, 0.965, barTitle, ...
       'Units','normalized', ...
       'HorizontalAlignment','center', ...
       'VerticalAlignment','top', ...
       'FontWeight','bold', ...
       'FontSize', 12, ...
       'Clipping','off');

  % --------- Daten: P + Phase ---------
  Q2 = squeeze(Q2);
  P  = real(0.5 * (Q2 .* conj(Q2)));

  A = P;
  A(~isfinite(A)) = -Inf;

  [vals, idx] = sort(A(:), "descend");
  [row, col]  = ind2sub(size(A), idx);

  good = isfinite(vals);
  vals = vals(good);
  row  = row(good);
  col  = col(good);

  K = min(barsPerRow, numel(vals));
  vals = vals(1:K);
  row  = row(1:K);
  col  = col(1:K);

  m0 = (size(Q2,1)+1)/2;
  labels = arrayfun(@(r,c) sprintf("M %d\nN %d", r-m0, c), row, col, ...
                    "UniformOutput", false);

  lin = sub2ind(size(Q2), row, col);
  ph  = angle(Q2(lin));

  cmap = viridis(256);
  nc = rows(cmap);

  vmin = min(vals);
  vmax = max(vals);
  if vmax <= vmin, vmax = vmin + 1; end

  ci = 1 + floor((vals - vmin) ./ (vmax - vmin) * (nc - 1));
  ci(ci<1)=1; ci(ci>nc)=nc;
  cols = cmap(ci,:);

  fsTick = 9;
  fsYLab = 11;

  draw_patchbars(axP,  vals, cols, labels, "P [W]",       fsTick, fsYLab);
  draw_patchbars(axPh, ph,   cols, labels, "Phase [rad]", fsTick, fsYLab);

end


function draw_patchbars(ax, y, cols, labels, ylab, fsTick, fsYLab)

  K = numel(y);
  cla(ax); hold(ax,'on');

  set(ax,'FontWeight','bold','FontSize',fsTick,'YGrid','on','LineWidth',1);
  grid(ax,'on'); box(ax,'on');

  ylabel(ax, ylab, 'FontWeight','bold','FontSize',fsYLab);

  xlim(ax,[0.5 K+0.5]);

  ymax = max(abs(y));
  if ymax==0, ymax=1; end
  if any(y<0)
    lim = ymax*1.25; ylim(ax,[-lim lim]);
  else
    ylim(ax,[0 ymax*1.25]);
  end

  w = 0.8;
  edgeCol = [0.25 0.25 0.30];

  for i = 1:K
    x0 = i - w/2;
    x1 = i + w/2;
    yi = y(i);

    patch('Parent',ax, ...
          'XData',[x0 x1 x1 x0], ...
          'YData',[0 0 yi yi], ...
          'FaceColor',cols(i,:), ...
          'EdgeColor',edgeCol, ...
          'LineWidth',0.35);
  end

  set(ax,'XTick',1:K,'XTickLabel',labels);
  set(ax,'XTickLabelRotation',90);

  hold(ax,'off');
end
