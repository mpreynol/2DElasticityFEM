%Set up Mesh Geometry:
[NN,NEL,X,Y] = GridRectangle(10,0.5,40,10);

%Set up Essential Boundary:
b1=[-eps,eps,-eps,2+eps,[0,0]];
b11=[10-eps,10+eps,-eps,2+eps,[0,0]];
BE=Boundary(NN,b1,b11);
[G,b]=Assemble.lagrange(BE);

% Set up Natural Boundary:
b2=[10-eps,10+eps,-eps,0.5+eps,[0,0]];
BN=Boundary(NN,b2); BN(BN==-Inf)=0;

% Set up Inputs:
Q=[0;-15];
C=Constit(100000,0.2,'Plane Strain').C;

% Set up Plotting Domain:
R=[0,15,-0.4,1];

%
%Set up Mesh Object as collection of element objects
Mesh=Element.empty(size(NEL,1),0);
for i=1:size(NEL,1)
    gNodes=NEL(i,:); x=NN(gNodes,2); y=NN(gNodes,3); 
    dof=reshape([NN(NEL(i,:),4),NN(NEL(i,:),5)]',[8,1]);
    h=BN(dof);
    Mesh(i)=Element(x,y,gNodes,dof,C,Q,h,2);
end
%
MeshPlot.plotOriginal(Mesh)
axis(R)

hold on
%
% Assembly Element and Force Vectors
[K,f]=Assemble.buildFromMesh(Mesh,size(NN,1)*2);

% % Solve System The Standard way:
% L=BE==-inf; % Indexes of unknown equations
% Kr=K(L,L); Br=BE(~L); fr=f(L); KRHS=K(L,~L); RHS=fr-KRHS*Br;
% ur=Kr\RHS;
% u=Assemble.reAssembleUnknowns(ur,BE);


% Solve System with Lagrange Multipliers:
[KA,fb]=Assemble.padLagrange(K,f,G,b);
ua=KA\fb;

%
% Populate solution back into Mesh Collection:
for i=1:size(Mesh,2)
    Mesh(i).u=ua(Mesh(i).dof);
    Mesh(i).setNodalResults();
end
%%
% Append Results for Node Array
NN=[NN,ua(NN(:,4)),ua(NN(:,5))];
for i=1:size(Mesh,2)
    NN(Mesh(i).nodes,8)=Mesh(i).E(1,:);
    NN(Mesh(i).nodes,9)=Mesh(i).E(2,:);
    NN(Mesh(i).nodes,10)=Mesh(i).E(3,:);
end

% Plot Deformed
MeshPlot.plotDeformed(Mesh,1)
axis(R)
hold on

%% Plot Results exx
Z=MeshPlot.buildSurface(X,Y,NN,8);
surface(X,Y,Z)
alpha(0.5)


%% Plot Results eyy
Z=MeshPlot.buildSurface(X,Y,NN,9);
surface(X,Y,Z)
alpha(0.5)


%% Plot Results exy
Z=MeshPlot.buildSurface(X,Y,NN,10);
surface(X,Y,Z)
alpha(0.5)





