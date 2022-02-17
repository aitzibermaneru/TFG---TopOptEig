classdef EulerBeamOptimizer < handle
    
    properties (Access = protected)
        optimizerType
    end
    
    properties (Access = private)
        nElem
        nConstraints
        length
        nValues
        youngModulus
        inertiaMoment
        minThick
        maxThick
        maxIter 
    end

    properties (Access = private)
        designVariable
        nIter
        mmaParams
    end
     
    methods (Access = public)
        
        function obj = EulerBeamOptimizer()
            obj.init()
            obj.computeVariablesMMA()
            obj.computeIterativeProcess()
        end

    end
    
    methods (Access = private)
        
        function init(obj)
            obj.nElem         = 10;
            obj.nConstraints  = 3; 
            obj.length        = 1/obj.nElem; 
            obj.nValues       = obj.nElem +1;
            obj.youngModulus  = 1;
            obj.inertiaMoment = 1;  
            obj.minThick      = 0.25;
            obj.maxThick      = 10;
            obj.optimizerType = 'MMA';
            obj.maxIter       = 1000;
        end

        function obj = computeVariablesMMA(obj)
            s.nElem         = obj.nElem;
            s.nConstraints  = obj.nConstraints; 
            s.length        = obj.length;
            s.youngModulus  = obj.youngModulus;
            s.inertiaMoment = obj.inertiaMoment;
            s.minThick      = obj.minThick;
            s.maxThick      = obj.maxThick;
            s.nValues       = obj.nValues;
            solution = MMAVariablesComputer(s);
            solution.compute();
            obj.designVariable = solution.designVariable;
            obj.mmaParams.xMin    = solution.xmin;
            obj.mmaParams.xMax    = solution.xmax;
            obj.mmaParams.xOld1   = solution.xold1;
            obj.mmaParams.xOld2   = solution.xold2;
            obj.nIter             = solution.loop;
            obj.mmaParams.uPP     = solution.upp;
            obj.mmaParams.lOW      = solution.low;
            obj.mmaParams.aMMA    = solution.a_mma;
            obj.mmaParams.a0Val   = solution.a0;
            obj.mmaParams.dVal    = solution.d;
            obj.mmaParams.cVal    = solution.c;
        end

        function obj = computeIterativeProcess(obj)
            s.mmaParams       = obj.mmaParams;
            s.nElem          = obj.nElem;
            s.nConstraints   = obj.nConstraints; 
            s.length         = obj.length;
            s.nValues        = obj.nValues;
            s.youngModulus   = obj.youngModulus;
            s.inertiaMoment  = obj.inertiaMoment;
            s.minThick       = obj.minThick;
            s.maxThick       = obj.maxThick;
            s.loop           = obj.nIter;
            s.designVariable = obj.designVariable;
            s.maxIter        = obj.maxIter;
            solution = IterativeProcessComputer(s);
            solution.compute();
            obj.designVariable = solution.designVariable;
        end     
        
    end
end