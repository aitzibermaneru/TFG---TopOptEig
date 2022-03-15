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

%% TO DO LIST
% 0. create eigenModes in Iterative; (DONE)
% 1. Cost to FirstEigenValueFuncitonal (shapeFunctional)   (DONE)
% 2. Cost is composed of FirstEigenValue  (DONE)
% 3. constraint to 3 shapeFunctionals (DONE)
% 4. constraint is composed of 3 shapeFunctionals ()
% 3. Use Optimizer_MMA with monitoring

% 4. plotBeam in 3D
%%
    
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
            obj.mmaParams = obj.createMMAparams(solution);
        end

        function cParams = createMMAparams(obj,s)
            cParams.xMin    = s.xmin;
            cParams.xMax    = s.xmax;
            cParams.xOld1   = s.xold1;
            cParams.xOld2   = s.xold2;
            obj.nIter       = s.loop;
            cParams.uPP     = s.upp;
            cParams.lOW     = s.low;
            cParams.aMMA    = s.a_mma;
            cParams.a0Val   = s.a0;
            cParams.dVal    = s.d;
            cParams.cVal    = s.c;
        end

        function obj = computeIterativeProcess(obj)
            s.mmaParams      = obj.mmaParams;
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
            s.optimizerType  = obj.optimizerType;
            solution = IterativeProcessComputer(s);
            solution.compute();
            obj.designVariable = solution.designVariable;
        end     
        
    end
end