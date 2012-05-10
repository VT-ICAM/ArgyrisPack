[xg yg wg] = gausspts();
qV = createArgyrisQuadratureValues(xg, yg);

someTriangle = [0.1190    0.9597    0.5853; 0.4984    0.3404    0.2238];

[C B b Th] = changeOfBasis(someTriangle);
% note that C2 is the transpose of C1 by the new convention.
[C2 B2 b2 Th2] = ArgyrisMaps(someTriangle);

[dx dy] = evalGradArgy(C, B, qV);
[dx2 dy2] = ArgyrisGradients(C2, B2, qV.derivativeX, qV.derivativeY);

[dxx dxy dyy] = evalHessArgy(C, Th, qV);
[dxx2 dxy2 dyy2] = ArgyrisHessians(C2, Th2, qV.derivativeXX, qV.derivativeXY, ...
                                   qV.derivativeYY);

disp('max difference in gradients:')
max([max(max(abs(dx2 - dx))), max(max(abs(dy2 - dy)))])
disp('max difference in second derivatives:')
max([max(max(dxx2 - dxx)), max(max(dxy2 - dxy)), max(max(dyy2 - dyy))])