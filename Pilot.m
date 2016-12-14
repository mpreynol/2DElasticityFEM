% These files are the modified Thermal Files for Elasticity
%Set up Mesh Geometry:
n=10;
Domain=[0,1,0,1]; %xmin,xmax,ymin,ymax
[NN,NEL,X,Y] = GridSquare(n,Domain(1),Domain(4)); % Set up Global Elements

%Set up Essential Boundary:
b1=[-eps,eps,Domain(3),Domain(4),[0,0]];
BE=Boundary(NN,b1);
%[G,b]=Assemble.lagrange(BE);

% Set up Natural Boundary:

% Set up Inputs:
Q=[0;-1];
C=Constit(100,0.2,'Plane Strain').C;

%%
%Set up Mesh Object as collection of element objects
Mesh=Element.empty(size(NEL,1),0);
for i=1:1
    dof=NEL(i,:); x=NN(dof,2); y=NN(dof,3); h=[zeros(1,4);zeros(1,4)]; 
    Mesh(i)=Element(x,y,dof,C,Q,h,2);
end
%%
% Assembly Element and Force Vectors
[K,f]=Assemble.buildFromMesh(Mesh,size(NN,1));

% % Solve System The Standard way:
% L=BE==-inf; % Indexes of unknown equations
% Kr=K(L,L); Br=BE(~L); fr=f(L); KRHS=K(L,~L); RHS=fr-KRHS*Br;
% ur=Kr\RHS;
% u=Assemble.reAssembleUnknowns(ur,BE);

% Solve System with Lagrange Multipliers:
[KA,fb]=Assemble.padLagrange(K,f,G,b);
ua=KA\fb;

%%
% Populate solution back into Mesh Collection:
for i=1:size(Mesh,2)
    Mesh(i).u=u(Mesh(i).dof);
end

% Append Results for Node Array for Parsing
NN=[NN,u];
Z=Assemble.buildSurface(X,Y,NN);

% Plot
figure(1)
surf(X,Y,Z,'facecolor','red','EdgeColor','red')
alpha(0.2)

    




