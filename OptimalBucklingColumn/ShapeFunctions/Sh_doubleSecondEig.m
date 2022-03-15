classdef Sh_doubleSecondEig < ShapeFunctional
    
    properties (Access = private)
        nElem
        designVariable
        stiffnessMatrix
        bendingMatrix
        elementalBendingMatrix
        D
        V
        v1
        v2
        lambda
    end

    properties (Access = private)
        eigModes
        bendingMat
        stiffnessMat
    end

    methods (Access = public)
        
        function obj = Sh_doubleSecondEig(cParams)
            obj.init(cParams)
        end

        function computeFunctionAndGradient(obj,iter)
            obj.computeFunction(iter);
            obj.computeGradient();
        end
        
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.nElem = cParams.nElem;
            obj.designVariable = cParams.designVariable;
            obj.eigModes = cParams.settings.eigMod;
            obj.bendingMat = cParams.settings.bendingMat;
            obj.stiffnessMat = cParams.settings.stiffnessMat;
        end
    end
    
    methods (Access = public)

        function computeFunction(obj,iter)
           [l] = obj.computeSettings(iter);
           x = obj.designVariable.value;
           N = obj.nElem;
           fx = x(N+1)-l(2);
           obj.value = fx;
        end

        function [l] = computeSettings(obj,iter)
            obj.stiffnessMat.compute();
            obj.stiffnessMatrix = obj.stiffnessMat.stiffnessMatrix;
            obj.bendingMat.compute();
            obj.elementalBendingMatrix = obj.bendingMat.elementalBendingMatrix;
            obj.bendingMatrix = obj.bendingMat.bendingMatrix;
            K = obj.stiffnessMatrix;
            B = obj.bendingMatrix;
            obj.eigModes.compute(K,B,iter);
            obj.lambda = obj.eigModes.lambda;
            obj.D      = obj.eigModes.D;
            obj.V      = obj.eigModes.V;
            obj.v1     = obj.eigModes.v1;
            obj.v2     = obj.eigModes.v2;
            obj.eigModes.compute(K,B,iter);
            l = obj.eigModes.lambda;
        end

        function computeGradient(obj,iter)
            Belem = obj.elementalBendingMatrix;
            x = obj.designVariable.value;
            N = obj.nElem;
            if abs(obj.D(2,2)-obj.D(1,1))> 1 % Simple eigenvalues
                W=zeros(2*N+2,2);
                for i=3:2*N
                    W(i,2)=obj.v2(i-2);
                end
                for i=1:N
                    dfdx(1,i)= -(2*x(i,1))*(W(2*(i-1)+1: 2*(i-1)+4,2)'*Belem*W(2*(i-1)+1: 2*(i-1)+4,2));
                end
            else % Dobles eigenvalues
                disp('dobles')
                obj.D
                Q1=zeros(2*N+2,1);
                Q2=zeros(2*N+2,1);
                dQ1=zeros(N,1);
                dQ2=zeros(N,1);
                    dQ1Q2=zeros(N,1);
                    for i=3:2*N
                        Q1(i,1)=obj.V(i-2,1);
                    end
                    for i=3:2*N
                        Q2(i,1)=obj.V(i-2,2);
                    end
                    for i=1:N
                        dQ1(i,1)= (2*x(i,1))*(Q1(2*(i-1)+1: 2*(i-1)+4,1)'*Belem*Q1(2*(i-1)+1: 2*(i-1)+4,1));
                        dQ2(i,1)= (2*x(i,1))*(Q2(2*(i-1)+1: 2*(i-1)+4,1)'*Belem*Q2(2*(i-1)+1: 2*(i-1)+4,1));
                        dQ1Q2(i,1)= (2*x(i,1))*(Q1(2*(i-1)+1: 2*(i-1)+4,1)'*Belem*Q2(2*(i-1)+1: 2*(i-1)+4,1));
                        A=[dQ1(i,1) dQ1Q2(i,1); dQ1Q2(i,1) dQ2(i,1)];
                        [U,R]=eigs(A,2,'SM');
                        S=sort(diag(R));
                        dfdx(1,i)=-S(2);
                    end
            end
            dfdx(1,N+1)=1;
            obj.gradient = dfdx';
        end
    end
    
end