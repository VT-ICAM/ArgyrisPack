function [dxx dxy dyy]=evalHessArgy(C, Th, quadratureValues)
%The following code was partially taken from Dominguez 2006

    dxx = quadratureValues.derivativeXX(:);
    dxy = quadratureValues.derivativeXY(:);
    dyy = quadratureValues.derivativeYY(:);
    k = length(quadratureValues.values);
    hess    = zeros(21,3*k);
    hess(:) = [dxx, dxy, dyy] * (inv(Th'));
    hess    = C' * hess;
    dxx=hess(:,1:k);
    dxy=hess(:,k+1:2*k);
    dyy=hess(:,2*k+1:3*k);

end
