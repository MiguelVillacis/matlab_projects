classdef TriMeshObj < handle
    %TRIMESH Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        g
        hmax
        Nodes
        Elements
    end
    
    methods
        function obj = TriMeshObj(g, hmax)
            %TRIMESH Construct an instance of this class
            obj.g = g;
            obj.hmax = hmax;
        end
        
        function obj = generateMesh(obj)
            mesh = gtrimesh(obj.g, obj.hmax);
            obj.Nodes = mesh.Nodes;
            obj.Elements = mesh.Elements;
        end
    end
end

