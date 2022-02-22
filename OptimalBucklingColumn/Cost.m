classdef Cost < handle
    
    properties (Access = public)
        gradient
        value
    end
    
    properties (Access = private)
        designVariable
        nElem
    end
    
    methods (Access = public)
        
        function obj = Cost(cParams)
            obj.init(cParams)
            obj.computeCost();
        end
        
        function computeFunctionAndGradient(obj)
            obj.computeCost();
        end

    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.designVariable = cParams.designVariable;
            obj.nElem          = cParams.nElem;
        end
        
        function computeCost(obj)
            N = obj.nElem;
            x = obj.designVariable.value;
            f0val = -x(N+1); 
            df0dx = zeros(N+1,1);
            df0dx(N+1) = -1;
            df0dx2 = 0*df0dx;
            obj.value = f0val;
            obj.gradient = df0dx;            
        end
        
    end
    
end