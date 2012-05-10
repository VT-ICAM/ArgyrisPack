function [x y w]=gausspts()
%determines the gauss quadrature points for the reference triangle
% this data came from John Burkardt 
% http://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html
%currently n is meaningless
% x is a row vector of x coordinates
% y is a row vector of y coordinates
% w is a column vector of weights

[p w] = simplexquad(5,2);

x = p(:,1)';
y = p(:,2)';
end
