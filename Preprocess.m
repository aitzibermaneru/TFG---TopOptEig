classdef Preprocess < handle

    properties (Access = public)
        E
        l
        I
        mmaParams
    end
    
    properties (Access = private)
       nel
       alpha
       beta
    end

    methods (Access = public)

        function obj = Preprocess (cParams)
                obj.init(cParams);
        end

        function obj = compute(obj)
            obj.computeConstants();
            obj.computeLimitations();
        end
    end

    methods (Access = private)

        function obj = init(obj,cParams)
                obj.nel   = cParams.nel;
                obj.alpha = cParams.alpha;
                obj.beta  = cParams.beta;
            end

        function computeConstants(obj)
            obj.E = 1;
            obj.l = 1/obj.nel;
            obj.I = 1;
        end

        function computeLimitations(obj)
            nElem  = obj.nel;
            alphav = obj.alpha;
            betav  = obj.beta;
            mma.n_val = nElem+1;
            n_v = mma.n_val;
            m = 3;
            x = ones(nElem+1,1);
            xminv = alphav*ones(nElem,1);
            mma.xmin = [xminv;0];
            xmaxv = betav*ones(nElem,1);
            mma.xmax = [xmaxv;0];
            mma.xold1=x;
            mma.xold2=x;
            mma.loop=0;
            mma.low = zeros(n_v,1);
            mma.upp = ones(n_v,1);
            mma.a0 = 1;
            mma.a_mma = zeros(m,1);
            mma.d = zeros(m,1);
            mma.c = 1000*ones(m,1);
            obj.mmaParams = mma;
        end
    end

end
