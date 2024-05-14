function y = Hill_line_neg(x,n,K,A)
% PIECEWISELINE   A line made of two pieces
% that is not continuous.

y = zeros(size(x));

y = A*1./((x.^n)+K);

end