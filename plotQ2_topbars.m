function h = plotQ2_topbars(Q2, varargin)
  % plotQ2_topbars(Q2, 'thr',0.001, 'maxBars',100, 'barsPerRow',25, 'y0',101)
  %


  % Defaults
  thr        = 0.001;
  maxBars    = 100;
  barsPerRow = 25;
  y0         = 101;

  % Style defaults
  titleFS = 16;   % Title font size
  axisFS  = 12;   % Axis label font size
  tickFS  = 10;   % Tick label font size
  valFS   = 9;    % Value text font size

  if mod(numel(varargin),2) ~= 0
    error("Optionen bitte als Key/Value Paare.");
  end
  for k = 1:2:numel(varargin)
    key = lower(varargin{k});
    val = varargin{k+1};
    switch key
      case {"thr","threshold"}
        thr = val;
      case {"maxbars","maxbar","max"}
        maxBars = val;
      case {"barsperrow","barspersrow","barsperline","perrow","perline"}
        barsPerRow = val;
      case {"y0","m0","mcenter","centerrow"}
        y0 = val;
      otherwise
        error("Unbekannte Option: %s", varargin{k});
    end
  end

  [m,n] = size(Q2);

  % --- find Peaks  ---
  mask = isfinite(Q2) & (Q2 >= thr);
  idx  = find(mask);

  if isempty(idx)
    fig = figure("Name","Q2 Top Bars","NumberTitle","off", "Color","w");
    axes("Parent", fig);
    text(0.1, 0.5, sprintf("Keine Werte >= %.4g", thr), "FontSize", 14, "FontWeight","bold");
    axis off;
    if nargout>0, h = struct("fig",fig); end
    return;
  end

  vals = Q2(idx);
  [valsSorted, ord] = sort(vals, "descend");
  idxSorted = idx(ord);

  K = min(maxBars, numel(idxSorted));
  idxSorted  = idxSorted(1:K);
  valsSorted = valsSorted(1:K);

  [ii, jj] = ind2sub([m,n], idxSorted);

  % --- Labels
  labels = cell(K,1);
  for t = 1:K
    mrel = ii(t) - y0;             % y0->0, y0-1->-1, ...
    labels{t} = sprintf("M%d,N%d", mrel, jj(t));
  end

  % --- Layout ---
  B = max(1, round(barsPerRow));
  rows = ceil(K / B);


  fig = figure("Name","Q2 Top Bars","NumberTitle","off", ...
               "Color","w", "MenuBar","figure", "ToolBar","figure");
  try
    set(fig, "Units","normalized", "Position",[0.05 0.08 0.90 0.84]);
  catch
  end


  try
    cmap = viridis(256);
  catch
    cmap = hot(256);
  end

  vminAll = min(valsSorted);
  vmaxAll = max(valsSorted);
  if vmaxAll <= vminAll, vmaxAll = vminAll + 1; end

  for r = 1:rows
    a = (r-1)*B + 1;
    b = min(r*B, K);

    ax = subplot(rows, 1, r, "Parent", fig);
    set(ax, "Color","w", "Box","off", ...
            "FontSize", tickFS, "FontWeight","bold", ...
            "XGrid","on","YGrid","on", ...
            "GridColor",[0.85 0.85 0.85], "GridAlpha", 1.0);

    x   = 1:(b-a+1);
    v   = valsSorted(a:b);
    lab = labels(a:b);


    bh = bar(ax, x, v);

    try
      set(bh, "FaceColor","flat", "EdgeColor","none");
      ci = 1 + floor((v - vminAll) ./ (vmaxAll - vminAll) * (rows(cmap)-1));
      ci = max(1, min(rows(cmap), ci));
      set(bh, "CData", cmap(ci,:));
    catch

      try, set(bh, "FaceColor",[0.2 0.4 0.7]); catch, end
    end

    % X-Ticks/Labels
    set(ax, "XTick", x);
    set(ax, "XTickLabel", lab);
    try
      set(ax, "XTickLabelRotation", 90);
    catch
    end


    ylabel(ax, "|Q2|", "FontSize", axisFS, "FontWeight","bold");


    vmax = max(v);
    if vmax <= 0, vmax = 1; end
    ylim(ax, [0, vmax * 1.28]);


    for t = 1:numel(v)
      text(ax, x(t), v(t) + 0.03*vmax, sprintf("%.2f", v(t)), ...
           "HorizontalAlignment","center", ...
           "VerticalAlignment","bottom", ...
           "FontSize", valFS, ...
           "FontWeight","bold", ...
           "Color",[0.15 0.15 0.15]);
    end


    if r == 1
      title(ax, sprintf("Top-%d Werte (thr ≥ %.4g)", K, thr), ...
            "FontSize", titleFS, "FontWeight","bold");
    end
  end


  han = axes("Position",[0 0 1 1], "Visible","off", "Parent", fig);
  xlabel(han, "Index (Mrel,N)", "FontSize", axisFS, "FontWeight","bold", "Visible","on");

  if nargout > 0
    h = struct("fig", fig);
  end
end
