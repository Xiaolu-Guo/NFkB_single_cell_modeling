function y = Hill_line_pos(x,n,K,A)
% PIECEWISELINE   A line made of two pieces
% that is not continuous.

y = zeros(size(x));

y = A*(x.^n)./((x.^n)+K);

end