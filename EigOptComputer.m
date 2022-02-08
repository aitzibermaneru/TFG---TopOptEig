classdef EigOptComputer < handle

        properties (Access = public)
            nel
            E
            I
            l
            mma
        end

        properties (Access = private)
          alpha 
          beta
        end

        methods (Access = public)

            function obj = EigOptComputer(cParams)
                obj.init(cParams);
            end

            function obj = compute(obj)
                obj.computePreprocess();
                obj.computeIterativeProcess();
            end

        end

        methods (Access = private)

            function obj = init(obj,cParams)
                obj.nel   = cParams.input.nel;
                obj.alpha = cParams.input.alpha;
                obj.beta  = cParams.input.beta;
            end

            function obj = computePreprocess(obj)
                s.nel   = obj.nel;
                s.alpha = obj.alpha;
                s.beta  = obj.beta;
                solution = Preprocess(s);
                solution.compute();
                obj.mma = solution.mmaParams;
                obj.l   = solution.l;
                obj.I   = solution.I;
                obj.E   = solution.E;
            end

            function obj = computeIterativeProcess(obj)
                s.nel   = obj.nel;
                s.E     = obj.E;
                s.I     = obj.I;
                s.l     = obj.l;
                s.alpha = obj.alpha;
                s.beta  = obj.beta;
                s.mma   = obj.mma;
                change = 1;
                while (change > 0.0005) & (s.mma.loop < 1000)
                    solution = IterativeProcessComputer(s);
                    solution.compute();
                end
            end

        end

end


