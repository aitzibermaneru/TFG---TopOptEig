classdef MatricesComputer < handle

    properties (Access = public)
        stiffnessMatrix
        bendingMatrix
    end

    properties (Access = private)
        E 
        l
        I 
        nel
        elementalStiffnessMatrix
        elementalBendingMatrix
    end

    methods (Access = public)
        function obj = MatricesComputer(cParams)
            obj.init(cParams);
        end

        function obj = compute(obj)
            obj.computeElementalMatrices();
            obj.computeGlobalMatrices();
        end
    end

    methods (Access = private)

        function obj = init(obj,cParams)
            obj.nel = cParams.nel;
            obj.E   = cParams.E;
            obj.I   = cParams.I;
            obj.l   = cParams.l;
        end

        function computeElementalMatrices(obj)
            lv = obj.l;
            Ev = obj.E;
            Iv = obj.I;
            Ke = obj.computeElementalStiffnessMatrix(lv);
            Be = obj.computeElementalBendingMatrix(lv,Ev,Iv);
            obj.elementalStiffnessMatrix = Ke;
            obj.elementalBendingMatrix   = Be;

        end

        function Ke = computeElementalStiffnessMatrix(obj,l)
            Ke = 1/(30*l)*[36 3*l -36 3*l; 
                3*l 4*l^2 -3*l -l^2; 
                -36 -3*l 36 -3*l; 
                3*l -l^2 -3*l 4*l^2];
        end

        function Be = computeElementalBendingMatrix(obj,l,E,I)
            Be = (E*I/(l^3))*[12 6*l -12 6*l;  
                6*l 4*l^2 -6*l 2*l^2; 
                -12 -6*l 12 -6*l; 
                6*l 2*l^2 -6*l 4*l^2]; 
        end

        function computeGlobalMatrices(obj)
            nElem = obj.nel;
            Ke = obj.elementalStiffnessMatrix;
            Be = obj.elementalBendingMatrix;
            B = sparse(2*nElem+2, 2*nElem+2);
            K = sparse(2*nElem+2, 2*nElem+2);
            x = ones(nElem+1,1);
            for iElem=1:nElem
                idof=[2*iElem-1; 2*iElem; 2*(iElem+1)-1; 2*(iElem+1)];
                B(idof,idof)=B(idof,idof)+(x(iElem)^2)*Be;
                K(idof,idof)=K(idof,idof)+ Ke;
            end
            obj.stiffnessMatrix = B;
            obj.bendingMatrix = K;
        end


    end
end