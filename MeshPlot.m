classdef MeshPlot<handle
    %MESHPLOT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Static)
        function plotOriginal(Mesh)
            for i=1:size(Mesh,2)
                x=Mesh(i).x; x=[x;x(1)];
                y=Mesh(i).y; y=[y;y(1)];
                lh.Color=[0,0,0,0.5]
                plot(x,y,'-',lh)
                hold on
            end
        end
        function plotDeformed(Mesh,scale)
             for i=1:size(Mesh,2)
                x=Mesh(i).x; x=[x;x(1)];
                dofx=Mesh(i).dof(1:2:8);
                dofy=Mesh(i).dof(2:2:8);
                ux=Mesh(i).u(1:2:8); ux=[ux;ux(1)];
                y=Mesh(i).y; y=[y;y(1)];
                uy=Mesh(i).u(2:2:8); uy=[uy;uy(1)];
                plot((x+scale*ux),(y+scale*uy),'-b')
                hold on
            end           
        end
    end
    
end

