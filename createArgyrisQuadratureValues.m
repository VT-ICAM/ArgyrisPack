% ==============================================================================
% File    : createArgyrisQuadratureValues.m
% Author  : David Wells < drwells at default Virginia Tech email >
% Purpose : wrapper around QuadratureValues to fill in the Argyris functions.
% ==============================================================================
function [quadValues] = createArgyrisQuadratureValues(xPoints, yPoints)
% -----------------------------------------------------------------------------
% function
%
% input  : xPoints, yPoints
%
% output : quadValues
%
% goal   : build QuadratureValues object.
% -----------------------------------------------------------------------------
    syms x y

    argyris = containers.Map('KeyType', 'double', 'ValueType', 'any');

    argyris(1) = 1 - 10*x^3 - 10*y^3 + 15*x^4 - 30*x^2*y^2 + 15*y^4 - 6*x^5 ...
        + 30*x^3*y^2 + 30*x^2*y^3 - 6*y^5;
    argyris(2) = 10*x^3 - 15*x^4 + 15*x^2*y^2 + 6*x^5 - 15*x^3*y^2 - 15*x^2*y^3;
    argyris(3) = 10*y^3 + 15*x^2*y^2 - 15*y^4 - 15*x^3*y^2 - 15*x^2*y^3 + 6*y^5;
    argyris(4) = x - 6*x^3 - 11*x*y^2 + 8*x^4 + 10*x^2*y^2 + 18*x*y^3 - 3*x^5 ...
        + x^3*y^2 - 10*x^2*y^3 - 8*x*y^4;
    argyris(5) = y - 11*x^2*y - 6*y^3 + 18*x^3*y + 10*x^2*y^2 + 8*y^4 - 8*x^4*y ...
        - 10*x^3*y^2 + x^2*y^3 - 3*y^5;
    argyris(6) = -4*x^3 + 7*x^4 - 3.5*x^2*y^2 - 3*x^5 + 3.5*x^3*y^2 + 3.5* ...
        x^2*y^3;
    argyris(7) = -5*x^2*y + 14*x^3*y + 18.5*x^2*y^2 - 8*x^4*y - 18.5*x^3*y^2 ...
        - 13.5*x^2*y^3;
    argyris(8) = -5*x*y^2 + 18.5*x^2*y^2 + 14*x*y^3 - 13.5*x^3*y^2 - 18.5*x^2*y^3 ...
        - 8*x*y^4;
    argyris(9) = -4*y^3 - 3.5*x^2*y^2 + 7*y^4 + 3.5*x^3*y^2 + 3.5*x^2*y^3 - ...
        3*y^5;
    argyris(10) = 0.5*x^2 - 1.5*x^3 + 1.5*x^4 - 1.5*x^2*y^2 - 0.5*x^5 + ...
        1.5*x^3*y^2 + x^2*y^3;
    argyris(11) = x*y - 4*x^2*y - 4*x*y^2 + 5*x^3*y + 10*x^2*y^2 + 5*x*y^3 - ...
        2*x^4*y - 6*x^3*y^2 - 6*x^2*y^3 - 2*x*y^4;
    argyris(12) = 0.5*y^2 - 1.5*y^3 - 1.5*x^2*y^2 + 1.5*y^4 + x^3*y^2 + ...
        1.5*x^2*y^3 - 0.5*y^5;
    argyris(13) = 0.5*x^3 - x^4 + 0.25*x^2*y^2 + 0.5*x^5 - 0.25*x^3*y^2 - ...
        0.25* x^2*y^3;
    argyris(14) = x^2*y - 3*x^3*y - 3.5*x^2*y^2 + 2*x^4*y + 3.5*x^3*y^2 + ...
        2.5*x^2*y^3;
    argyris(15) = 1.25*x^2*y^2 - 0.75*x^3*y^2 - 1.25*x^2*y^3;
    argyris(16) = 1.25*x^2*y^2 - 1.25*x^3*y^2 - 0.75*x^2*y^3;
    argyris(17) = x*y^2 - 3.5*x^2*y^2 - 3*x*y^3 + 2.5*x^3*y^2 + 3.5*x^2*y^3 ...
        + 2*x*y^4;
    argyris(18) = 0.5*y^3 + 0.25*x^2*y^2 - y^4 - 0.25*x^3*y^2 - 0.25*x^2*y^3 ...
        + 0.5*y^5;
    argyris(19) = 16*x^2*y - 32*x^3*y - 32*x^2*y^2 + 16*x^4*y + 32*x^3*y^2 + ...
        16*x^2*y^3;
    argyris(20) = -16*x*y^2 + 32*x^2*y^2 + 32*x*y^3 - 16*x^3*y^2 - 32*x^2*y^3 ...
        - 16*x*y^4;
    argyris(21) = sqrt(2) * (8*x^2*y^2 - 8*x^3*y^2 - 8*x^2*y^3);

    quadValues = QuadratureValues(argyris, xPoints, yPoints);
end