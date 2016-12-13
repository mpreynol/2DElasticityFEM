function [B] = Boundary(NN,varargin)
% Boundary Function takes in Nodal Array and a collection of boundary
% values builds an array NN elements long specifying the assigned
% boundary value or -inf if DOF is free.

% Assumptions
%a. Boundary is the closed set
%b. NN is formatted as per "GridSquare" function
%c. To specify an edge specify a min and max correspoding a thin threshold


% Syntax for Boundary Specification:
% [xmin, xmax, ymin,ymax,value]
B=ones(size(NN,1),1)*-inf;
for i=1:nargin-1
    infoArray=varargin{i};
    for j=1:length(NN)
        if (NN(j,2)>=infoArray(1) && NN(j,2)<=infoArray(2) && NN(j,3)>=infoArray(3) && NN(j,3)<=infoArray(4))
            B(j)=infoArray(5);
        end
    end
end
end

