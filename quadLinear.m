classdef quadLinear < handle
    %Class defines an object containing shape functions for a linear quad
    %element
    % The intent is that this object is created per element and computes
    % shape function at the gauss integration points
	% See Page 148 of "The Finite Element Method - Thomas J.R Hughes"
	% Written by Mathew Reynolds Nov 10, 2016
    
    properties
        x=[]; % Nodal x cordinates of the element
        y=[]; % Nodal y cordinates of the element
        X
        Y
        N=[]; % The Shape function array that is used for the forcing arrays
        NE=[] % Shape Functions compiled for 2D Elasticity
        Nxi=[];
        Neta=[];
        Xxi;
        Xeta;
        Yxi;
        Yeta;
        j;
        Nx=[];
        Ny=[];
        B=[]; % The gradient array that is used for the stiffness matrix
        BE=[]; % Shape Function Gradient for 2D Elasticity
    end
    
    methods
        % Overloaded Constructor:
        function obj = quadLinear(x,y)
            obj.x=x;
            obj.y=y;
            obj.N=zeros(1,4);
            obj.NE=zeros(2,8);
            obj.Nxi=zeros(1,4);
            obj.Neta=zeros(1,4);
            obj.Nx=zeros(1,4);
            obj.Ny=zeros(1,4);
            obj.B=zeros(2,4);
            obj.BE=zeros(3,8);
        end
        % Set Method to reset coordinates 
        function setCords(obj,x,y)
            obj.x=x;
            obj.y=y;
        end
        
        % Set All Properties based on xi and eta
        function setAll(obj,xi,eta)
            obj.N(1)=0.25*(1-xi)*(1-eta);
            obj.N(2)=0.25*(1+xi)*(1-eta);
            obj.N(3)=0.25*(1+xi)*(1+eta);
            obj.N(4)=0.25*(1-xi)*(1+eta);
            obj.NE=[obj.N(1)*eye(2),obj.N(2)*eye(2),obj.N(3)*eye(2),obj.N(4)*eye(2)];
            obj.Nxi(1)=0.25*(eta-1);
            obj.Nxi(2)=0.25*(1-eta);
            obj.Nxi(3)=0.25*(1+eta);
            obj.Nxi(4)=-0.25*(1+eta);
            obj.Neta(1)=0.25*(xi-1);
            obj.Neta(2)=-0.25*(1+xi);
            obj.Neta(3)=0.25*(1+xi);
            obj.Neta(4)=0.25*(1-xi);
            obj.setX();
            obj.setY();
            obj.setXxi();
            obj.setXeta();
            obj.setYxi();
            obj.setYeta();
            obj.setJ();
            obj.setNx();
            obj.setNy();
            obj.setB();
            obj.setBE();
        end
    
        % Set Method for x
        function setX(obj)
            obj.X=0;
            for i=1:length(obj.N)
                obj.X=obj.X+obj.N(i)*obj.x(i);
            end
        end
        % Set Method for y
        function setY(obj)
            obj.Y=0;
            for i=1:length(obj.N)
                obj.Y=obj.Y+obj.N(i)*obj.y(i);
            end
        end
        % Set Method for xXi
        function setXxi(obj)
            obj.Xxi=0;
            for i=1:length(obj.N)
                obj.Xxi=obj.Xxi+obj.Nxi(i)*obj.x(i);
            end
        end      
        % Set Method for xEta
        function setXeta(obj)
            obj.Xeta=0;
            for i=1:length(obj.N)
                obj.Xeta=obj.Xeta+obj.Neta(i)*obj.x(i);
            end
        end
        % Set Method for yXi
        function setYxi(obj)
            obj.Yxi=0;
            for i=1:length(obj.N)
                obj.Yxi=obj.Yxi+obj.Nxi(i)*obj.y(i);
            end
        end      
        % Set Method for yEta
        function setYeta(obj)
            obj.Yeta=0;
            for i=1:length(obj.N)
                obj.Yeta=obj.Yeta+obj.Neta(i)*obj.y(i);
            end 
        end    
        
        % Set Method for j
        function setJ(obj)
            obj.j=obj.Xxi*obj.Yeta-obj.Xeta*obj.Yxi;
        end
        
        % Set Method for Nx
        function setNx(obj)
            for i=1:length(obj.N)
                obj.Nx(i)=(obj.Nxi(i)*obj.Yeta-obj.Neta(i)*obj.Yxi)/obj.j;
            end
        end
        
        % Set Method for Ny
        function setNy(obj)
             for i=1:length(obj.N)
                obj.Ny(i)=-(obj.Nxi(i)*obj.Xeta-obj.Neta(i)*obj.Xxi)/obj.j;
            end           
        end

        % Set Method for B
        function setB(obj)
            obj.B=[obj.Nx;obj.Ny];
        end
        
        % Set Method for BE
        function setBE(obj)
            obj.BE(1,:)=[obj.Nx(1),0,obj.Nx(2),0,obj.Nx(3),0,obj.Nx(4),0];
            obj.BE(2,:)=[0,obj.Ny(1),0,obj.Ny(2),0,obj.Ny(3),0,obj.Ny(4)];
            obj.BE(3,:)=[obj.Ny(1),obj.Nx(1),obj.Ny(2),obj.Nx(2),obj.Ny(3),obj.Nx(3),obj.Ny(4),obj.Nx(4)];
        end
        
        %Get Method for xi,eta
        function [R]=getXiEta(obj,X,Y)
            % Xi and Eta are solved for by using Newton Raphson
            % See Page 53 of Numerical Renaissance for derivation of
            % multivariate NR for nonlinear equations
            % Inital Guess:
            R=[0.1;0.1];
            tol = 0.0000001;
            h=[1;1];
            %NewtonRaphson Until we solve for xi,eta
            while norm(h)>=tol
               obj.setAll(R(1),R(2));
               A=-[obj.x';obj.y']*obj.B';
               f=[X;Y]-[obj.N*obj.x;obj.N*obj.y];
               h=A\-f;
               R=R+h;
            end
        end
    end
    
    
end

