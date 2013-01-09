argyris = containers.Map('KeyType', 'double', 'ValueType', 'any');
argyrisdx = containers.Map('KeyType', 'double', 'ValueType', 'any');
argyrisdy = containers.Map('KeyType', 'double', 'ValueType', 'any');
argyrisdxx = containers.Map('KeyType', 'double', 'ValueType', 'any');
argyrisdxy = containers.Map('KeyType', 'double', 'ValueType', 'any');
argyrisdyy = containers.Map('KeyType', 'double', 'ValueType', 'any');

syms x y

% functions for the triangle xs = [1/2, 3/2, 1], ys = [1, 1, 3/2]
argyris(1) = -6*x^5 - 60*x^3*y^2 - 120*x^2*y^3 - 30*x*y^4 + 30*x^4 + 24*y^5 + 120*x^3*y + 600*x^2*y^2 + 360*x*y^3 - 115*x^3 - 150*y^4 - 840*x^2*y - 1185*x*y^2 + 405*x^2 + 290*y^3 + 1410*x*y - 120*y^2 - 4575/8*x - 150*y + 875/8;
argyris(2) = 6*x^5 + 60*x^3*y^2 - 120*x^2*y^3 + 30*x*y^4 - 30*x^4 + 24*y^5 - 120*x^3*y + 240*x^2*y^2 + 120*x*y^3 + 115*x^3 - 210*y^4 - 120*x^2*y - 495*x*y^2 - 45*x^2 + 530*y^3 + 510*x*y - 570*y^2 - 1185/8*x + 270*y - 371/8;
argyris(3) = 240*x^2*y^3 - 48*y^5 - 840*x^2*y^2 - 480*x*y^3 + 360*y^4 + 960*x^2*y + 1680*x*y^2 - 360*x^2 - 820*y^3 - 1920*x*y + 690*y^2 + 720*x - 120*y - 62;
argyris(4) = -3*x^5 - 28*x^3*y^2 - 46*x^2*y^3 + 7*x*y^4 + 31/2*x^4 + 22*y^5 + 56*x^3*y + 248*x^2*y^2 + 50*x*y^3 - 115/2*x^3 - 315/2*y^4 - 358*x^2*y - 321*x*y^2 + 723/4*x^2 + 809/2*y^3 + 464*x*y - 448*y^2 - 3335/16*x + 405/2*y - 725/32;
argyris(5) = -8*x^4*y - 2*x^3*y^2 + 2*x^2*y^3 - 22*x*y^4 + 8*x^4 - 18*y^5 + 38*x^3*y - 9*x^2*y^2 + 100*x*y^3 - 36*x^3 + 139*y^4 - 38*x^2*y - 311/2*x*y^2 + 45*x^2 - 807/2*y^3 + 255/2*x*y + 2205/4*y^2 - 50*x - 725/2*y + 375/4;
argyris(6) = -3*x^5 - 28*x^3*y^2 + 46*x^2*y^3 + 7*x*y^4 + 29/2*x^4 - 22*y^5 + 56*x^3*y - 80*x^2*y^2 - 134*x*y^3 - 107/2*x^3 + 287/2*y^4 + 22*x^2*y + 335*x*y^2 + 129/4*x^2 - 641/2*y^3 - 296*x*y + 322*y^2 + 1289/16*x - 293/2*y + 785/32;
argyris(7) = -8*x^4*y + 2*x^3*y^2 + 2*x^2*y^3 + 22*x*y^4 + 8*x^4 - 18*y^5 + 26*x^3*y - 21*x^2*y^2 - 108*x*y^3 - 28*x^3 + 95*y^4 - 2*x^2*y + 431/2*x*y^2 + 21*x^2 - 391/2*y^3 - 351/2*x*y + 753/4*y^2 + 46*x - 167/2*y + 55/4;
argyris(8) = -68*x^3*y^2 + 20*x*y^4 + 136*x^3*y + 204*x^2*y^2 - 116*x*y^3 - 68*x^3 - 20*y^4 - 408*x^2*y + 41*x*y^2 + 204*x^2 + 116*y^3 + 186*x*y - 177*y^2 - 131*x + 86*y - 5;
argyris(9) = -28*x^2*y^3 - 20*y^5 + 98*x^2*y^2 + 56*x*y^3 + 114*y^4 - 112*x^2*y - 196*x*y^2 + 42*x^2 - 279*y^3 + 224*x*y + 727/2*y^2 - 84*x - 246*y + 135/2;
argyris(10) = -1/2*x^5 - 5/2*x^3*y^2 - 9/2*x^2*y^3 + x*y^4 + 11/4*x^4 + 5/2*y^5 + 5*x^3*y + 95/4*x^2*y^2 + 7/2*x*y^3 - 33/4*x^3 - 18*y^4 - 34*x^2*y - 231/8*x*y^2 + 163/8*x^2 + 375/8*y^3 + 173/4*x*y - 859/16*y^2 - 685/32*x + 105/4*y - 225/64;
argyris(11) = -2*x^4*y - 3*x^3*y^2 - x^2*y^3 - x*y^4 + 2*x^4 - y^5 + 15*x^3*y + 23/2*x^2*y^2 + 6*x*y^3 - 12*x^3 + 15/2*y^4 - 69/2*x^2*y - 77/4*x*y^2 + 24*x^2 - 87/4*y^3 + 137/4*x*y + 257/8*y^2 - 20*x - 105/4*y + 75/8;
argyris(12) = 1/2*x^3*y^2 - 1/2*x^2*y^3 - 5/2*x*y^4 - 3/2*y^5 - x^3*y + 1/4*x^2*y^2 + 27/2*x*y^3 + 1/2*x^3 + 49/4*y^4 + x^2*y - 201/8*x*y^2 - 3/4*x^2 - 305/8*y^3 + 79/4*x*y + 903/16*y^2 - 45/8*x - 40*y + 175/16;
argyris(13) = 1/2*x^5 + 5/2*x^3*y^2 - 9/2*x^2*y^3 - x*y^4 - 9/4*x^4 + 5/2*y^5 - 5*x^3*y + 35/4*x^2*y^2 + 29/2*x*y^3 + 25/4*x^3 - 16*y^4 - 4*x^2*y - 289/8*x*y^2 - 25/8*x^2 + 287/8*y^3 + 131/4*x*y - 583/16*y^2 - 291/32*x + 67/4*y - 181/64;
argyris(14) = 2*x^4*y - 3*x^3*y^2 + x^2*y^3 - x*y^4 - 2*x^4 + y^5 - x^3*y + 13/2*x^2*y^2 + 2*x*y^3 + 4*x^3 - 11/2*y^4 - 15/2*x^2*y - 37/4*x*y^2 + 55/4*y^3 + 49/4*x*y - 125/8*y^2 - 4*x + 31/4*y - 11/8;
argyris(15) = -1/2*x^3*y^2 - 1/2*x^2*y^3 + 5/2*x*y^4 - 3/2*y^5 + x^3*y + 13/4*x^2*y^2 - 23/2*x*y^3 - 1/2*x^3 + 29/4*y^4 - 5*x^2*y + 145/8*x*y^2 + 9/4*x^2 - 105/8*y^3 - 47/4*x*y + 179/16*y^2 + 21/8*x - 9/2*y + 11/16;
argyris(16) = 9*x^2*y^3 - 5*y^5 - 59/2*x^2*y^2 - 18*x*y^3 + 65/2*y^4 + 32*x^2*y + 59*x*y^2 - 23/2*x^2 - 299/4*y^3 - 64*x*y + 619/8*y^2 + 23*x - 71/2*y + 43/8;
argyris(17) = 6*x^3*y^2 + 2*x*y^4 - 12*x^3*y - 18*x^2*y^2 - 6*x*y^3 + 6*x^3 - 2*y^4 + 36*x^2*y + 45/2*x*y^2 - 18*x^2 + 6*y^3 - 35*x*y - 21/2*y^2 + 33/2*x + 11*y - 9/2;
argyris(18) = x^2*y^3 + 3*y^5 - 7/2*x^2*y^2 - 2*x*y^3 - 35/2*y^4 + 4*x^2*y + 7*x*y^2 - 3/2*x^2 + 165/4*y^3 - 8*x*y - 393/8*y^2 + 3*x + 59/2*y - 57/8;
argyris(19) = 16*x^4*y - 32*x^2*y^3 - 16*x^4 + 16*y^5 - 64*x^3*y + 128*x^2*y^2 + 64*x*y^3 + 64*x^3 - 112*y^4 - 72*x^2*y - 256*x*y^2 - 24*x^2 + 280*y^3 + 272*x*y - 304*y^2 - 80*x + 145*y - 25;
argyris(20) = -64*sqrt(1/2)*x^3*y^2 - 64*sqrt(1/2)*x^2*y^3 + 64*sqrt(1/2)*x*y^4 + 64*sqrt(1/2)*y^5 + 128*sqrt(1/2)*x^3*y + 416*sqrt(1/2)*x^2*y^2 - 192*sqrt(1/2)*x*y^3 - 480*sqrt(1/2)*y^4 - 64*sqrt(1/2)*x^3 - 640*sqrt(1/2)*x^2*y - 48*sqrt(1/2)*x*y^2 + 1328*sqrt(1/2)*y^3 + 288*sqrt(1/2)*x^2 + 416*sqrt(1/2)*x*y - 1672*sqrt(1/2)*y^2 - 240*sqrt(1/2)*x + 960*sqrt(1/2)*y - 200*sqrt(1/2);
argyris(21) = -64*sqrt(1/2)*x^3*y^2 + 64*sqrt(1/2)*x^2*y^3 + 64*sqrt(1/2)*x*y^4 ...
    - 64*sqrt(1/2)*y^5 + 128*sqrt(1/2)*x^3*y - 32*sqrt(1/2)*x^2*y^2 - 448*sqrt(1/2)*x*y^3 + 352*sqrt(1/2)*y^4 - 64*sqrt(1/2)*x^3 - 128*sqrt(1/2)*x^2*y + 848*sqrt(1/2)*x*y^2 - 688*sqrt(1/2)*y^3 + 96*sqrt(1/2)*x^2 - 608*sqrt(1/2)*x*y + 616*sqrt(1/2)*y^2 + 144*sqrt(1/2)*x - 256*sqrt(1/2)*y + 40*sqrt(1/2);

for i=1:21
    argyrisdx(i) = diff(argyris(i), x);
    argyrisdy(i) = diff(argyris(i), y);
    argyrisdxx(i) = diff(argyris(i), x, 2);
    argyrisdxy(i) = diff(diff(argyris(i), y), x);
    argyrisdyy(i) = diff(argyris(i), y, 2);
end

% generate comparison points.
xs = rand(50,1);
ys = (1 - xs).*rand(50,1);

[C B b Th] = apGlobalMapsMex([1/2, 3/2, 1], [1, 1, 3/2]);

functionRef            = apLocalFunctionsMex(xs, ys);
[dxRef dyRef]          = apLocalGradientsMex(xs, ys);
[dxxRef dxyRef dyyRef] = apLocalHessiansMex(xs, ys);

values = apGlobalFunctionsMex(C, functionRef);
[dx dy] = apGlobalGradientsMex(C, B, dxRef, dyRef);
[dxx dxy dyy] = apGlobalHessiansMex(C, Th, dxxRef, dxyRef, dyyRef);

% map the comparison points to global coordinates.
globalPoints = B*[xs';ys'];
globalxs = globalPoints(1,:) + b(1);
globalys = globalPoints(2,:) + b(2);

valuesSymbolic = zeros(size(values));
dxSymbolic = zeros(size(dx));
dySymbolic = zeros(size(dy));
dxxSymbolic = zeros(size(dxx));
dxySymbolic = zeros(size(dxy));
dyySymbolic = zeros(size(dyy));

for i=1:21
    for j=1:size(valuesSymbolic, 2)
        valuesSymbolic(i,j) = subs(subs(argyris(i), x, globalxs(j)), y, globalys(j));
        dxSymbolic(i,j)     = subs(subs(argyrisdx(i), x, globalxs(j)), y, globalys(j));
        dySymbolic(i,j)     = subs(subs(argyrisdy(i), x, globalxs(j)), y, globalys(j));
        dxxSymbolic(i,j)    = subs(subs(argyrisdxx(i), x, globalxs(j)), y, globalys(j));
        dxySymbolic(i,j)    = subs(subs(argyrisdxy(i), x, globalxs(j)), y, globalys(j));
        dyySymbolic(i,j)    = subs(subs(argyrisdyy(i), x, globalxs(j)), y, globalys(j));
    end
end

disp('Error Results:')
disp(['Max Error in function values : ' max(max(abs(values - ...
                                                  valuesSymbolic)))]);
disp(['Max Error in dx              : ' max(max(abs(dx - dxSymbolic)))]);
disp(['Max Error in dy              : ' max(max(abs(dy - dySymbolic)))]);
disp(['Max Error in dxx             : ' max(max(abs(dxx - dxxSymbolic)))]);
disp(['Max Error in dxy             : ' max(max(abs(dxy - dxySymbolic)))]);
disp(['Max Error in dyy             : ' max(max(abs(dyy - dyySymbolic)))]);