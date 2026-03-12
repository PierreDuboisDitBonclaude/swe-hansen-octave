function update_plots_cb2(h, Q)


  s_lo = strtrim(get(h.e_q_lo,'String'));
  s_hi = strtrim(get(h.e_q_hi,'String'));
  if isempty(s_lo) || isempty(s_hi), return; end

  q_lo = str2double(s_lo);
  q_hi = str2double(s_hi);
  if any(isnan([q_lo q_hi])), return; end
  if q_lo >= q_hi, return; end


  mode_chart4(Q,1,3,q_lo,q_hi,'pow',h.axPlot1,'Q_{1mn} TE-Modes');
  mode_chart4(Q,2,3,q_lo,q_hi,'pow',h.axPlot2,'Q_{2mn} TM-Modes');
  bars4(h.pBars1,squeeze(Q(1,:,:)),20,'Q_{1mn} TE-Modes');
  bars4(h.pBars2,squeeze(Q(2,:,:)),20,'Q_{2mn} TM-Modes');

  R = qsmn_stats(Q);


  Nshow = min(200, size(R,1));
  R = R(1:Nshow,:);


  data = cell(Nshow,5);
  sw=0;
  for i=1:Nshow

    if (R(i,7)>= 99.995 && sw==0 )
        data{i,1} = ['*' sprintf('%d|%d|%d', R(i,1), R(i,2), R(i,3)) ];
        data{i,2} = ['*' sprintf('%.2f', R(i,4)) ];
        data{i,3} = ['*' sprintf('%.2f', R(i,5)) ];
        data{i,4} = ['*' sprintf('%.2f', R(i,6)) ];
        data{i,5} = ['*' sprintf('%.2f',R(i,7)) ];
        sw=1;
    elseif (R(i,7)>= 99.995 || sw==0 )
        data{i,1} = sprintf('%d|%d|%d', R(i,1), R(i,2), R(i,3));
        data{i,2} = sprintf('%.2f', R(i,4));
        data{i,3} = sprintf('%.2f', R(i,5));
        data{i,4} = sprintf('%.2f', R(i,6));
        data{i,5} = sprintf('%.2f',R(i,7));
    end
  end

  set(h.tbl,'Data',data);

end
