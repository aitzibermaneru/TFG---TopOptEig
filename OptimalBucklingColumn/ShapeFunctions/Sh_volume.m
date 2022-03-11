classdef Sh_volume < ShapeFunctional
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = Sh_volume(cParams)
            obj.init(cParams)

        end
        function computeFunctionAndGradient(obj)
            obj.computeFunction();
            obj.computeGradient();
        end

    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            
        end

        function computeFunction(obj)
            fx = (1/N)*sum(x(1:N))-1;
            obj.value = fx;
        end

        function computeGradient(obj)
            dfdx(3,1:N)=(1/N)*ones(1,N);
            obj.gradient = dfdx;
        end

    end

end