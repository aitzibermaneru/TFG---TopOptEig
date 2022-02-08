classdef IterativeProcessComputer < handle 

    
    properties (Access = private)
        E
        l
        I
        mmaParams
        nel
        alpha
        beta
        stiffnessMatrix
        bendingMatrix
    end

    methods (Access = public)
        function obj = IterativeProcessComputer(cParams)
            obj.init(cParams);
        end

        function obj = compute(obj)
            obj.computeFEM();
            obj.computeBoundary();
            obj.computeEig();
            obj.computeMMA();
            obj.plotting();
        end
    end

    methods (Access = private)

        function obj = init(obj,cParams)
            obj.nel       = cParams.nel;
            obj.E         = cParams.E;
            obj.I         = cParams.I;
            obj.l         = cParams.l;
            obj.alpha     = cParams.alpha;
            obj.beta      = cParams.beta;
            obj.mmaParams = cParams.mma;
        end

        function obj = computeFEM(obj)
            s.nel = obj.nel;
            s.l   = obj.l;
            s.E   = obj.E;
            s.I   = obj.I;
            solution = MatricesComputer(s);
            solution.compute();
            obj.stiffnessMatrix = solution.stiffnessMatrix;
            obj.bendingMatrix   = solution.bendingMatrix;
        end

        function computeBoundary(obj)

        end
        
        function computeEig(obj)
            
        end

        function computeMMA(obj)

        end

        function plotting(obj)

        end

    end
end