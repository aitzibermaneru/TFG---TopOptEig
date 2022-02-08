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
            [c1,c2,c3,c4,c5] = obj.coeffsStiffness(l);
            Ke(1,1:4) = c1*[c2 c3 -c2 c3];
            Ke(2,1:4) = c1*[c3 c4 -c3 -c5];
            Ke(3,1:4) = c1*[-c2 -c3 c2 -c3];
            Ke(4,1:4) = c1*[c3 -c5 -c3 c4];
%           Ke = 1/(30*l)*[36 3*l -36 3*l; 
%                 3*l 4*l^2 -3*l -l^2; 
%                 -36 -3*l 36 -3*l; 
%                 3*l -l^2 -3*l 4*l^2];
        end

        function [c1,c2,c3,c4,c5] = coeffsStiffness(obj,l)
            c1 = 1/(30*l);
            c2 = 36;
            c3 = 3*l;
            c4 = 4*l^2;
            c5 = l^2;
        end

        function Be = computeElementalBendingMatrix(obj,l,E,I)
            [c1,c2,c3,c4] = obj.coeffsBending(l,E,I);
            Be(1,1:4) = c1*[c2 c3 -c2 c3];
            Be(2,1:4) = c1*[c3 c4 -c3 c4/2];
            Be(3,1:4) = c1*[-c2 -c3 c2 -c3];
            Be(4,1:4) = c1*[c3 c4/2 -c3 c4];
%             Be = (E*I/(l^3))*[12 6*l -12 6*l;  
%                 6*l 4*l^2 -6*l 2*l^2; 
%                 -12 -6*l 12 -6*l; 
%                 6*l 2*l^2 -6*l 4*l^2]; 
        end

        function [c1,c2,c3,c4] = coeffsBending(obj,l,E,I)
            c1 = E*I/l^3;
            c2 = 12;
            c3 = 6*l;
            c4 = 4*l^2;
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
            obj.stiffnessMatrix = K;
            obj.bendingMatrix = B;
        end

    end
end