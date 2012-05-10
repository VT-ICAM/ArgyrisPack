function [C,B,b,Th] = changeOfBasis(p)
% This function calculates the required transformations for the change of basis
% from khat to k.
% p - 3 x 2 vector containing the vertices of the triangle K
% C - Transformation for the Argyris Basis Function from K to Khat
% B - Transformation for the gradient in K to Khat
% b - point corresponding to (x1,y1)
% Th - Transformation for the hessian in K to Khat

b = p(:,1);
%Affine Transformations
B = [p(:,2)-p(:,1), p(:,3)-p(:,1)]; 
Th = [B(1,1)^2, 2*B(1,1)*B(2,1), B(2,1)^2;
  B(1,2)*B(1,1), B(1,2)*B(2,1)+B(1,1)*B(2,2), B(2,1)*B(2,2);
  B(1,2)^2, 2*B(2,2)*B(1,2), B(2,2)^2];

%vectors
v = [p(:,2)-p(:,1), p(:,3)-p(:,1), p(:,3)-p(:,2)]; 
%side lengths
L = diag([norm(v(:,1)), norm(v(:,2)), norm(v(:,3))]);
A = L^(-2)*[ 0 1; -1 0; -1/sqrt(2) -1/sqrt(2)]*B';
%Rotation Matrix, rotates counter clockwise pi/2
R = [0 -1; 1 0];

%Now calculate the matrix D
f = sum(A' .* (R*v));
g = sum(A' .* v); % this is equivalent to dot(A', v) but is much faster.
D = zeros(21,24);

% manual copying is unfortunately far faster than blkdiag.
D(1:3,1:3) = eye(3);
D(4:5,4:5) = B';
D(6:7,6:7) = B';
D(8:9,8:9) = B';
D(10:12,10:12) = Th;
D(13:15,13:15) = Th;
D(16:18,16:18) = Th;
D(19:21,19:21) = diag(f);
D(19:21,22:24) = diag(g);

% Here [diag(f) diag(g)] corresponds to the matrix Q

%Now calculate the Matrix E
w = [v(1,:).^2; 2*v(1,:).*v(2,:); v(2,:).^2]';
E = zeros(24,21);

E(1:18,1:18) = eye(18);
E(19:21, 19:21) = L;
E(22:24,1:3) = 15/8*[-1 1 0; -1 0 1; 0 -1 1];
E(22:24,4:9) = -7/16*[v(1,1), v(2,1), v(1,1), v(2,1), 0, 0; 
                      v(1,2), v(2,2), 0,      0,      v(1,2), v(2,2); 
                      0,      0,      v(1,3), v(2,3), v(1,3), v(2,3)];
E(22:24,10:18) = 1/32*[-w(1,1), -w(1,2), -w(1,3), w(1,1), w(1,2), w(1,3), 0, 0, 0;
                       -w(2,1), -w(2,2), -w(2,3), 0, 0, 0, w(2,1), w(2,2), w(2,3);
                        0, 0, 0, -w(3,1), -w(3,2), -w(3,3), w(3,1), w(3,2) w(3,3)];

C=D*E;

end
