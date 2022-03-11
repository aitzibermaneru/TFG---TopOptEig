classdef Cost < CC
    
    properties (Access = private)
        designVariable
        nElem
    end
    
    methods (Access = public)
        
        function obj = Cost(cParams)
            obj.init(cParams)
        end
        
        function computeFunctionAndGradient(obj)
            obj.computeFunctions();
            obj.computeGradients();
        end

    end
    
end