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
    end

    properties (Access = private)
        designVariable
        xMin
        xMax
        xOld1
        xOld2
        looping
        lOW
        uPP
        a0Val
        aMMA
        dVal
        cVal
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
            obj.xMin    = solution.xmin;
            obj.xMax    = solution.xmax;
            obj.xOld1   = solution.xold1;
            obj.xOld2   = solution.xold2;
            obj.looping = solution.loop;
            obj.uPP     = solution.upp;
            obj.lOW     = solution.low;
            obj.aMMA    = solution.a_mma;
            obj.a0Val   = solution.a0;
            obj.dVal    = solution.d;
            obj.cVal    = solution.c;
        end

        function obj = computeIterativeProcess(obj)
            s.nElem          = obj.nElem;
            s.nConstraints   = obj.nConstraints; 
            s.length         = obj.length;
            s.nValues        = obj.nValues;
            s.youngModulus   = obj.youngModulus;
            s.inertiaMoment  = obj.inertiaMoment;
            s.minThick       = obj.minThick;
            s.maxThick       = obj.maxThick;
            s.loop           = obj.looping;
            s.designVariable = obj.designVariable;
            s.xMin           = obj.xMin;
            s.xMax           = obj.xMax;
            s.xOld1          = obj.xOld1;
            s.xOld2          = obj.xOld2;
            s.lOW            = obj.lOW;
            s.uPP            = obj.uPP;
            s.a0Val          = obj.a0Val;
            s.aMMA           = obj.aMMA;
            s.dVal           = obj.dVal;
            s.cVal           = obj.cVal;
            solution = IterativeProcessComputer(s);
            solution.compute();
            obj.designVariable = solution.designVariable;
        end     
        
    end
end