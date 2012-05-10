% ==============================================================================
% File    : QuadratureValues.m
% Author  : David Wells < drwells at default Virginia Tech email >
% Purpose : Class for precomputing quadrature values of basis functions and
%           various derivatives.
% ==============================================================================

% ------------------------------------------------------------------------------
% Summary of input arguments
% ------------------------------------------------------------------------------
% functions - containers.Map() of symbolic objects, indexed by integers
%             starting at 1.
% x         - x-values of quadrature points
% y         - y-values of quadrature points

classdef QuadratureValues
    properties (Access = public)
        values
        derivativeX
        derivativeY
        derivativeXX
        derivativeXY
        derivativeYY
    end
    % --------------------------------------------------------------------------
    % Public Methods
    % --------------------------------------------------------------------------
    methods (Access = public)
        function self = QuadratureValues(functions, xPoints, yPoints)
            % constructor
            % force points to be single columns.
            xPoints = xPoints(:);
            yPoints = yPoints(:);

            syms x y

            functionXDerivatives  = containers.Map('KeyType', 'double', ...
                                                        'ValueType', 'any');
            functionYDerivatives  = containers.Map('KeyType', 'double', ...
                                                        'ValueType', 'any');
            functionXXDerivatives = containers.Map('KeyType', 'double', ...
                                                        'ValueType', 'any');
            functionXYDerivatives = containers.Map('KeyType', 'double', ...
                                                        'ValueType', 'any');
            functionYYDerivatives = containers.Map('KeyType', 'double', ...
                                                        'ValueType', 'any');

            % calculate all the derivatives as functions.
            for i=1:length(functions)
                functionXDerivatives(i)  = ...
                    matlabFunction(diff(functions(i),x));
                functionYDerivatives(i)  = ...
                    matlabFunction(diff(functions(i),y));
                functionXXDerivatives(i) = ...
                    matlabFunction(diff(functions(i),x,2));
                functionXYDerivatives(i) = ...
                    matlabFunction(diff(diff(functions(i),x),y));
                functionYYDerivatives(i) = ...
                    matlabFunction(diff(functions(i),y,2));
            end

            % convert the given functions from symbolic statements to
            % anonymous functions.
            for i=1:length(functions)
                functions(i) = matlabFunction(functions(i));
            end

            % create the lookup tables. Rows correspond to functions and
            % columns correspond to quadrature points.
            self.values       = zeros(length(functions), length(xPoints));
            self.derivativeX  = zeros(length(functions), length(xPoints));
            self.derivativeY  = zeros(length(functions), length(xPoints));
            self.derivativeXX = zeros(length(functions), length(xPoints));
            self.derivativeXY = zeros(length(functions), length(xPoints));
            self.derivativeYY = zeros(length(functions), length(xPoints));

            for i=1:length(functions)
                self.values(i,:)      = feval(functions(i), ...
                                              xPoints, yPoints);
                self.derivativeX(i,:) = feval(functionXDerivatives(i), ...
                                              xPoints, yPoints);
                self.derivativeY(i,:) = feval(functionYDerivatives(i), ...
                                              xPoints, yPoints);
                self.derivativeXX(i,:) = feval(functionXXDerivatives(i), ...
                                               xPoints, yPoints);
                self.derivativeXY(i,:) = feval(functionXYDerivatives(i), ...
                                               xPoints, yPoints);
                self.derivativeYY(i,:) = feval(functionYYDerivatives(i), ...
                                               xPoints, yPoints);
            end
        end
    end % end methods
end % end classdef
