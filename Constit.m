classdef Constit<handle
    %CONSTIT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        lambda; % Lames Constant
        lambdaOver; % Modified Lames Constand for Plane Stress
        Mu; % Lames Constant
        nu; % Poissons Ratio
        E; % Youngs Modulus
        alpha; % Constiutive Switching property
        token; % Plane Strain or Plane Stress Not case or space sensitive
        C; % Constituitive Tensor
    end
    
    methods
        % Overloaded Constructor
        function obj = Constit(E,nu,token)
            obj.E=E;
            obj.nu=nu;
            obj.token=lower(strtrim(token));
            obj.lambda=nu*E/((1+nu)*(1-2*nu));
            obj.Mu=E/(2*(1+nu));
            obj.lambdaOver=obj.lambda*2*obj.Mu/(obj.lambda+2*obj.Mu);
            obj.C=zeros(3,3);
            obj.setTensor()
        end
        
        %Method that returns 3x3 constiutive tensor
        function setTensor(obj)
           
            switch obj.token
                case 'planestrain'
                     obj.C(1,:)=[obj.lambda + 2*obj.Mu, obj.lambda, 0];
                     obj.C(2,:)=[obj.lambda, obj.lambda + 2*obj.Mu, 0];
                     obj.C(3,:) = [0,0,obj.Mu];
                case 'planestress'
                     obj.C(1,:)=[obj.lambdaOver + 2*obj.Mu, obj.lambdaOver, 0];
                     obj.C(2,:)=[obj.lambdaOver, obj.lambdaOver + 2*obj.Mu, 0];
                     obj.C(3,:) = [0,0,obj.Mu];
            end
        end
        
    end
    
end

