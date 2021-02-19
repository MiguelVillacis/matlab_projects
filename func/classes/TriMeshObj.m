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
            mmesh = gtrimesh(obj.g, obj.hmax);
            obj.Nodes = mmesh.Mesh.Nodes;
            obj.Elements = mmesh.Mesh.Elements;
        end
    end
end

