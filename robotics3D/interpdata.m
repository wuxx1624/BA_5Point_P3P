function [y_new] = interpdata(y_old, x_old, x_new)
% this function interpolates data y_old, which is sampled at times x_old,
% to obtain the values y_new at times x_new using simple linear
% interpolation


y_new = zeros(length(x_new), size(y_old,2));

parfor i = 1:length(x_new)
  
  tprev = find(x_old < x_new(i), 1, 'last');
  tnext = find(x_old >= x_new(i), 1, 'first');
  
  if isempty(tprev)
    y_new(i,:) = y_old(tnext,:);
  elseif isempty(tnext)
    y_new(i,:) = y_old(tprev,:);
  else
    dx = x_old(tnext) - x_old(tprev);
    dy = y_old(tnext,:) - y_old(tprev,:);
    
    dtau = x_new(i) - x_old(tprev);
    
    y_new(i,:) = y_old(tprev,:) + dtau * dy/dx;
  end
  
  
end