classdef CC < handle
    
    properties (Access = public)
        value
        gradient
    end
    
    properties (GetAccess = public, SetAccess = private)
        nSF
        shapeFunctions
    end
    
    properties (Access = private)
        nElem
        sizeDesigVar
        designVariable
    end
    
    methods (Access = protected)
        
        function obj = init(obj,cParams)
            obj.nSF = 0;
            obj.nElem = cParams.nElem;
            obj.sizeDesigVar = size(cParams.designVariable.value); 
            obj.createShapeFunctions(cParams);
            obj.designVariable = cParams.designVariable;
        end
        
    end
    
    methods (Access = protected)

        function createShapeFunctions(obj,cParams)
            nS = cParams.nShapeFunction; % cost 1, constraint 3
            for iS = 1:nS
                s.designVariable = cParams.designVariable; 
                s.type = cParams.type{iS};
                s.nElem = obj.nElem;
                shapeFunction = ShapeFunctional.create(s);
                obj.append(shapeFunction);
            end
        end
    end

    methods (Access = public)

        function computeFunctions(obj)
        %    obj.initValueAndGradient();
            for iSF = 1:length(obj.shapeFunctions)
                obj.shapeFunctions{iSF}.computeFunction();
                obj.value(iSF) = obj.shapeFunctions{iSF}.value;
            end
        end    

        function computeGradients(obj)
            for iSF = 1:length(obj.shapeFunctions)
                obj.shapeFunctions{iSF}.computeGradient();
                obj.gradient(:,iSF) = obj.shapeFunctions{iSF}.gradient;
            end
        end

        function append(obj,shapeFunction)
            obj.shapeFunctions{obj.nSF+1} = shapeFunction;
            obj.nSF = obj.nSF+1;
        end
    end

    methods (Access = private)

        function initValueAndGradient(obj)
            obj.value = 0;
            obj.gradient = zeros(obj.sizeDesigVar);  % no coge bien el size          
        end

    end
    
end