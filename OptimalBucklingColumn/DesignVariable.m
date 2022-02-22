classdef DesignVariable < handle
    
    properties (Access = public)
       value
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = DesignVariable()
            obj.init()            
        end
        
        function update(obj,x)
            obj.value = x;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            
        end
        
    end
    
end