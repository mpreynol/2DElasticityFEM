%Set up Mesh Geometry:
[NN,NEL,X,Y] = GridRectangle(20,1,20,2);

%Set up Essential Boundary:
b1=[-eps,eps,-eps,+eps,[0,-Inf]];
b2=[-eps,eps,0.5-eps,0.5+eps,[0,0]];
b3=[-eps,eps,1-eps,1+eps,[0,-Inf]];
BE=Boundary(NN,b1,b2,b3);
%[G,b]=Assemble.lagrange(BE);

% Set up Natural Boundary:
b2=[20-eps,20+eps,-eps,1+eps,[0,1]];
BN=Boundary(NN,b2); BN(BN==-Inf)=0;

% Set up Inputs:
Q=[0;0];
C=Constit(1E6,0.2,'Plane Stress').C;

% Set up Plotting Domain:
R=[0,22,-2,2];

%
%Set up Mesh Object as collection of element objects
Mesh=Element.empty(size(NEL,1),0);
for i=1:size(NEL,1)
    gNodes=NEL(i,:); x=NN(gNodes,2); y=NN(gNodes,3); 
    dof=reshape([NN(NEL(i,:),4),NN(NEL(i,:),5)]',[8,1]);
    h=BN(dof);
    Mesh(i)=Element(x,y,gNodes,dof,C,Q,h,2,@parabolicStress);
end
%
MeshPlot.plotOriginal(Mesh)
axis(R)

hold on
%
% Assembly Element and Force Vectors
[K,f]=Assemble.buildFromMesh(Mesh,size(NN,1)*2);

% Solve System The Standard way:
L=BE==-inf; % Indexes of unknown equations
Kr=K(L,L); Br=BE(~L); fr=f(L); KRHS=K(L,~L); RHS=fr-KRHS*Br;
ur=Kr\RHS;
u=Assemble.reAssembleUnknowns(ur,BE);


% Solve System with Lagrange Multipliers:
% [KA,fb]=Assemble.padLagrange(K,f,G,b);
% ua=KA\fb;

%
% Populate solution back into Mesh Collection:
for i=1:size(Mesh,2)
    Mesh(i).u=u(Mesh(i).dof);
    Mesh(i).setNodalResults();
end
%
% Append Results for Node Array
NN=[NN,u(NN(:,4)),u(NN(:,5))];
for i=1:size(Mesh,2)
    NN(Mesh(i).nodes,8)=Mesh(i).sigma(1,:);
    NN(Mesh(i).nodes,9)=Mesh(i).sigma(2,:);
    NN(Mesh(i).nodes,10)=Mesh(i).sigma(3,:);
end

% Plot Deformed
MeshPlot.plotDeformed(Mesh,1)
axis(R)
hold on

%% Plot Results sxx
Z=MeshPlot.buildSurface(X,Y,NN,8);
surface(X,Y,Z)
alpha(0.5)



%% Plot Results syy
Z=MeshPlot.buildSurface(X,Y,NN,9);
surface(X,Y,Z)
alpha(0.5)


%% Plot Results sxy
Z=MeshPlot.buildSurface(X,Y,NN,10);
surface(X,Y,Z)
alpha(0.5)

%% Section Cut
ySample=0:0.01:1;
xSpot=10;
% Assemble Data Arrays
for w=1:length(ySample)
    for o=1:size(Mesh,2)
        if (sum(Mesh(o).y>=(ySample(w)-eps))>0 && sum(Mesh(o).y<=(ySample(w))+eps)>0 &&...
                Mesh(o).x(1)==xSpot)
            uSample(w)=Mesh(o).getU(xSpot,ySample(w));
        end
    end
end

%% Plot
for c=1:length(ySample)
   value=uSample(c).dsigma(3);
   plot(ySample(c),value,'.');
   hold on
end




