function [ NN,NEL,X,Y ] = GridSquare(n,xmin,xmax)

% Function creates a square finite Element Grid with total width = w and
% n^2 elements
dw=(xmax-xmin)/n; % Node Spacing
nN=(n+1)^2; % Number of Nodes
[X,Y]=meshgrid(xmin:dw:xmax);

% Define Arrays:
NN= [(1:nN)',reshape(Y,size(Y,1)^2,1),reshape(X,size(X,1)^2,1)]; % Global Nodes
NEL=zeros(n^2,4);

for i=1:n^2 % Loops through the elements
   row=ceil(i/(n));
   column=mod(i-1,(n))+1;
   NEL(i,:)=[row*(n+1)+column-(n+1),row*(n+1)+column-(n+1)+1,(row+1)*(n+1)+column-(n+1)+1,(row+1)*(n+1)+column-(n+1)];
end

end

