classdef Assemble < handle
    % Assembly Class consists of static methods to assemble Stiffness
    % Matrix and Force Vectors from a mesh collection
    
    properties
    end
    
    methods(Static)
        function [K,f]= buildFromMesh(Mesh,n)
            % Method Loops through a mesh of size n and returns a stiffness
            % and force vector
           K=zeros(n,n);
           f=zeros(n,1);
           for kitten=1:length(Mesh)
                dof=Mesh(kitten).dof; % Use dof array as logical indexer
                K(dof,dof)=K(dof,dof)+Mesh(kitten).K;
                f(dof)=f(dof)+Mesh(kitten).f;    
           end
        end
        function [ufull] = reAssembleUnknowns(ureduced,BE)
            % Method reassembles a full 'u' vector for a reduced oned
            L=BE==-inf;
            ufull=zeros(length(BE),1);
            counter=1;
            for i=1:length(BE)
               if L(i)==1
                   ufull(i)=ureduced(counter);
                   counter=counter+1;
               else
                   ufull(i)=BE(i);
               end
            end
        end
        function [Z]=buildSurface(X,Y,NN)
            % Method Builds a surface from X and Y meshgrids a a return
            % nodal array
            Z=zeros(size(X,1),size(X,2));
            for i =1: length(X)
                for j=1: length(Y)
                    for k=1:length(NN)
                        if (NN(k,2)==X(i,j) && NN(k,3)==Y(i,j))
                            Z(i,j)=NN(k,4);
                        end
                    end
                end
            end
        end
        
        function [G,b]=lagrange(BE)
            % Method Builds Lagrange Multiplier Matrices and vectors 
            % Input is a array of length "NN" which lists the specified
            % nodal value or -Inf
            index=1:size(BE); % Logical indexer
            liveNodes=index(BE~=-Inf);
            G=zeros(length(liveNodes),length(BE));
            b=zeros(length(liveNodes),1);
            for i =1:length(liveNodes)
               G(i,liveNodes(i))=1;
               b(i)=BE(liveNodes(i));
            end
        end
        
        function [KA,fb]=padLagrange(K,f,G,b)
           % Function pads lagrange matrix and vector onto our 'k'
           n=size(K,1); nA = size(G,1);
           assert(n==size(G,2))
           assert(nA==length(b))
           KA=[K,G';G,zeros(nA,nA)];
           fb=[f;b];
        end
    end
end
    


