function [dx dy]=evalGradArgy(C,B, quadratureValues)
%The following code was partially taken from Dominguez 2006

    k=length(quadratureValues.values);
    zx = quadratureValues.derivativeX(:);
    zy = quadratureValues.derivativeY(:);
    grads=zeros(21,2*k);
    grads(:) = [zx zy] * (inv(B));
    grads=C'*grads; %now transform to K and evaluate the actual gradiant
    dx=grads(:,1:k);
    dy=grads(:,k+1:2*k);
end
