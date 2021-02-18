classdef Result
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        t_coord;
        Nodof;
        NNdx;
        NNdy;
        Dof;
        Loadx;
        LoadNdx;
        Nel;
        Node_Total;
        NNtotal;
        scalex;
        scaley;
        KK;
        TablaQ;
        Fnodal;
        TablaDef;
        TableRes;
    end
    
    methods
        function obj = untitled3(inputArg1,inputArg2)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

