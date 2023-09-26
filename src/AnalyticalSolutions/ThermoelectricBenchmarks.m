classdef ThermoelectricBenchmarks
    %THERMOELECTRICBENCHMARKS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        l1
        l2
        seebeck
        k
        sigma
        area
    end
    
    methods
        function obj = ThermoelectricBenchmarks(reader,mesh)
            %THERMOELECTRICBENCHMARKS Construct an instance of this class
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

