%Set up Mesh Geometry:
[NN,NEL,X,Y] = GridRectangle(10,0.5,50,5);

%Set up Essential Boundary:
b1=[-eps,eps,-eps,2+eps,[0,0]];
b11=[10-eps,10+eps,-eps,2+eps,[0,0]];
BE=Boundary(NN,b1,b11);
[G,b]=Assemble.lagrange(BE);

% Set up Natural Boundary:
b2=[10-eps,10+eps,-eps,0.5+eps,[0,0]];
BN=Boundary(NN,b2); BN(BN==-Inf)=0;

% Set up Inputs:
Q=[0;-10];
C=Constit(100000,0.2,'Plane Strain').C;

% Set up Plotting Domain:
R=[0,12,-1,1];

%
%Set up Mesh Object as collection of element objects
Mesh=Element.empty(size(NEL,1),0);
for i=1:size(NEL,1)
    gNodes=NEL(i,:); x=NN(gNodes,2); y=NN(gNodes,3); 
    dof=reshape([NN(NEL(i,:),4),NN(NEL(i,:),5)]',[8,1]);
    h=BN(dof);
    Mesh(i)=Element(x,y,dof,C,Q,h,2);
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
end

% Append Results for Node Array
NN=[NN,u(NN(:,4)),u(NN(:,5))];

% Plot Deformed
MeshPlot.plotDeformed(Mesh,1)
axis(R)
hold on
    




