classdef Element < handle
    %ELEMENT2DElasticity assembles LHS and RHS vectors by gauss int
    % Element is currently quad linear elements
	% Written by Mathew Reynolds Dec 12, 2016
    
    properties
        K=[]; % Stiffness Matrix
        dof=[]; % List of DOF for record keeping
        x=[]; % List of x cordinates for element
        y=[]; % List of y cordinates for element
        orderInt; % Order of integration to develop components
        G2=[]; % Guassian Integration Array in 2D
        G1=[]; % Guassian Integration Array in 1D
        Shape; % object containing shape functions
        fQ=[]; % Body Forcing Term
        fh=[]; % Flux Forcing Term
        f=[]; % Total Forcing Term
        C=[]; %  Constituitive Tensor
        Q=[]; % Body Force Term
        u=[]; % Deformations at the nodes once solved
        L2=0; % Value of L2 norm for the element
        h; % input boundary data.
    end
    
    methods
        % Overloaded constructor with nodal inputs and integration order
        function obj = Element(x,y,dof,C,Q,h,orderInt)
            obj.x=x;
            obj.y=y;
            obj.dof=dof;
            obj.C=C;
            obj.Q=Q;
            obj.h=h;
            obj.orderInt=orderInt;
            obj.G2=Quadrature.twoDim(orderInt);
            obj.G1=Quadrature.oneDim(orderInt);
            obj.Shape=quadLinear(x,y); % Initialize Shape Function
            obj.K=zeros(8,8); % Initalize Stiffess
            obj.setStiffness();
            obj.fQ=zeros(8,1); % Initalize Body Force
            obj.setBodyForce();
            obj.fh=zeros(8,1); % Initalize Flux Force
            obj.setTraction();
            obj.f=obj.fQ+obj.fh;
        end
        
        function setStiffness(obj)
            for i=1:size(obj.G2,1) % perform Guass Integration [-1,1] over domain
                obj.Shape.setAll(obj.G2(i,1),obj.G2(i,2));
                obj.K=obj.K+obj.Shape.BE'*obj.C*obj.Shape.BE*obj.G2(i,3)*obj.Shape.j;
            end
        end
        
        function setBodyForce(obj)
            % Need to rewrite forcing functions to use mapping
            for i=1:size(obj.G2,1) % perform Guass Integration [-1,1] over domain
                obj.Shape.setAll(obj.G2(i,1),obj.G2(i,2));
                obj.fQ=obj.fQ+obj.Shape.NE'*obj.Q*obj.Shape.j*obj.G2(i,3);             
            end            
        end
        
        function setTraction(obj)
            if sum((sum(obj.h~=0))) % Then we have tractions
                if obj.h(2,1)~=0 && obj.h(2,2)~=0 % Then this edge has tractions
                    n=[0;-1];
                    % Surface 1: eta=-1
                    for i=1:size(obj.G1,1) % perform Guass Integration [-1,1] over domain
                        obj.Shape.setAll(obj.G1(i,1),-1);
                        js=sqrt((obj.Shape.Xxi)^2+(obj.Shape.Yxi)^2);
                        obj.fh=obj.fh+obj.Shape.N'*dot(obj.h*obj.Shape.N',n')*js*obj.G1(i,2);
                    end
                end
                if obj.h(1,2)~=0 && obj.h(1,3)~=0 % Then this edge has tractions
                    n=[1;0];
                    % Surface 2: xi=1
                    for i=1:size(obj.G1,1) % perform Guass Integration [-1,1] over domain
                        obj.Shape.setAll(1,obj.G1(i,1));
                        js=sqrt((obj.Shape.Xeta)^2+(obj.Shape.Yeta)^2);
                        obj.fh=obj.fh+obj.Shape.N'*dot(obj.h*obj.Shape.N',n')*js*obj.G1(i,2);
                    end
                end
                if obj.h(2,3)~=0 && obj.h(2,4)~=0 % Then this edge has tractions
                    n=[0;1];
                    % Surface 3: eta=1
                    for i=1:size(obj.G1,1) % perform Guass Integration [-1,1] over domain
                        obj.Shape.setAll(obj.G1(i,1),1);
                        js=sqrt((obj.Shape.Xxi)^2+(obj.Shape.Yxi)^2);
                        obj.fh=obj.fh+obj.Shape.N'*dot(obj.h*obj.Shape.N',n')*js*obj.G1(i,2);
                    end
                end
                if obj.h(1,4)~=0 && obj.h(1,1)~=0 % Then this edge has tractions
                    n=[-1;0];
                    % Surface 4: xi=-1
                    for i=1:size(obj.G1,1) % perform Guass Integration [-1,1] over domain
                        obj.Shape.setAll(-1,obj.G1(i,1));
                        js=sqrt((obj.Shape.Xeta)^2+(obj.Shape.Yeta)^2);
                        obj.fh=obj.fh+obj.Shape.N'*dot(obj.h*obj.Shape.N',n')*js*obj.G1(i,2);
                    end
                end
            end
    end
        
        function setL2(obj)
           obj.L2=0;
           L2G2=Quadrature.twoDim(5);
           for i=1:size(L2G2,1)
               obj.Shape.setAll(L2G2(i,1),L2G2(i,2));
               uh=obj.Shape.N*obj.u;
               uexact=sinsin(obj.Shape.X,obj.Shape.Y);
               obj.L2=obj.L2+((uh-uexact)^2)*L2G2(i,3)*obj.Shape.j;
           end
           
        end
        
        function interU=getU(obj,X,Y)
            % Method returns the interpolated u at spatial coordinates X,Y
            [natural]=obj.Shape.getXiEta(X,Y);
            obj.Shape.setAll(natural(1),natural(2));
            interU=obj.Shape.N*obj.u;    
        end
    end
    
end

