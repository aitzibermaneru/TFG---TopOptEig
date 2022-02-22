classdef MMAVariablesComputer < handle
    
    properties (Access = public)
        designVariable
        xmin
        xmax
        xold1
        xold2
        loop
        low
        upp
        a0
        a_mma
        d
        c
    end
    
    properties (Access = private) 
        nElem
        nConstraints
        nValues
        length
        youngModulus
        inertiaMoment
        minThick
        maxThick
    end
    
    methods (Access = public)
        
        function obj = MMAVariablesComputer(cParams)
            obj.init(cParams)
        end
        
        function obj = compute(obj)
            obj.computeVariablesMMA();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.nElem        = cParams.nElem;
            obj.nConstraints = cParams.nConstraints; % m
            obj.length       = cParams.length;
            obj.youngModulus = cParams.youngModulus;
            obj.inertiaMoment = cParams.inertiaMoment;
            obj.minThick      = cParams.minThick;
            obj.maxThick      = cParams.maxThick;
            obj.nValues      =  cParams.nValues;
        end
        
        function computeVariablesMMA(obj)
            n_val = obj.nValues;
            m = obj.nConstraints;
            alpha = obj.minThick;
            beta = obj.maxThick;
            N = obj.nElem;
            x=ones(N+1,1);
            obj.xmin=alpha*ones(N,1);
            obj.xmin=[obj.xmin; 0];
            obj.xmax=beta*ones(N,1);
            obj.xmax=[obj.xmax; 1000];
            obj.xold1=x;
            obj.xold2=x;
            obj.loop=0;
            obj.low = zeros(n_val,1);
            obj.upp = ones(n_val,1);
            obj.a0 = 1;
            obj.a_mma = zeros(m,1);
            obj.d = zeros(m,1);
            obj.c = 1000*ones(m,1);
            obj.designVariable = x;
        end
        
    end
    
end