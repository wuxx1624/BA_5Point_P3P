% pseudo-hubert
function [d_mark] = phubert(d, b)

C = 2*b^2 * ( sqrt(1 + (abs(d) / b).^2) - 1);

w = sqrt(C) ./ abs(d);

d_mark = w .* d;

end