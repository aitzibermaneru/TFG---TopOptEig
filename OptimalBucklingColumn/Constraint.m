classdef Constraint < handle
    
    properties (GetAccess = public, SetAccess = private)
        mode1
        mode2 
        D
    end
    
    properties (Access = public)
        gradient
        value
    end
    
    properties (Access = private)
        bendingMatrix
        stiffnessMatrix
        v1
        v2
        V
        lambda
        Be
    end
    
    properties (Access = private)
       designVariable 
       nElem
       freeNodes
       length
       youngModulus
       inertiaMoment
       nConstraints
    end
    
    methods (Access = public)
        
        function obj = Constraint(cParams)
            obj.init(cParams)
            obj.computeElementalBendingMatrix();
            obj.computeStiffnessMatrix();            
        end
        
        function computeFunctionAndGradient(obj)
            obj.computeBendingMatrix();
            obj.computeEigenModesAndValues();
            obj.computeBucklingModes();            
            obj.computeConstraintFunction();
            obj.computeConstraintDerivative();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.designVariable = cParams.designVariable;
            obj.nElem = cParams.nElem;
            obj.freeNodes = cParams.freeNodes;
            obj.length = cParams.length;
            obj.youngModulus = cParams.youngModulus;
            obj.inertiaMoment = cParams.inertiaMoment;   
            obj.nConstraints = cParams.nConstraints;
        end
        
        function computeBendingMatrix(obj)
            x = obj.designVariable.value;
            N = obj.nElem;
            B=sparse(2*N+2, 2*N+2);    
            for iElem = 1: N
                iDof=[2*iElem-1; 2*iElem; 2*(iElem+1)-1; 2*(iElem+1)];
                B(iDof,iDof)=B(iDof,iDof)+(x(iElem)^2)*obj.Be;
            end
            obj.bendingMatrix = B;
        end        
        
        function computeEigenModesAndValues(obj)            
            free = obj.freeNodes;
            Bfree = obj.bendingMatrix(free,free);
            Kfree = obj.stiffnessMatrix(free,free);
            [V,D] = obj.computeEigenFunctionAndValues(Bfree,Kfree);
            lambda=sort(diag(D));
            [v1,v2] = obj.reorderModes(lambda,V,D);
            obj.v1 = v1;
            obj.v2 = v2;
            obj.V  = V;
            obj.D  = D;
            obj.lambda = lambda;
        end    
        
        function K = computeStiffnessMatrix(obj)
            Ke = obj.computeElementalStiffnessMatrix();            
            N = obj.nElem;
            K = sparse(2*N+2, 2*N+2);     
            for iElem = 1: N
                iDof=[2*iElem-1; 2*iElem; 2*(iElem+1)-1; 2*(iElem+1)];
                K(iDof,iDof)=K(iDof,iDof)+ Ke;
            end
            obj.stiffnessMatrix = K;
        end  
        
        function Be = computeElementalBendingMatrix(obj)
            L = obj.length;
            E = obj.youngModulus;
            I = obj.inertiaMoment;
            [c1,c2,c3,c4] = obj.coeffsBending(L,E,I);
            Be(1,1:4) = c1*[c2 c3 -c2 c3];
            Be(2,1:4) = c1*[c3 c4 -c3 c4/2];
            Be(3,1:4) = c1*[-c2 -c3 c2 -c3];
            Be(4,1:4) = c1*[c3 c4/2 -c3 c4];
            obj.Be = Be;
        end

        function [c1,c2,c3,c4] = coeffsBending(obj,L,E,I)
            c1 = E*I/L^3;
            c2 = 12;
            c3 = 6*L;
            c4 = 4*L^2;
        end

        function Ke = computeElementalStiffnessMatrix(obj)
            L = obj.length;       
            [c1,c2,c3,c4,c5] = obj.coeffsStiffness(L);
            Ke(1,1:4) = c1*[c2 c3 -c2 c3];
            Ke(2,1:4) = c1*[c3 c4 -c3 -c5];
            Ke(3,1:4) = c1*[-c2 -c3 c2 -c3];
            Ke(4,1:4) = c1*[c3 -c5 -c3 c4];
        end

        function [c1,c2,c3,c4,c5] = coeffsStiffness(obj,L)
            c1 = 1/(30*L);
            c2 = 36;
            c3 = 3*L;
            c4 = 4*L^2;
            c5 = L^2;
        end        
        
        function fx = computeConstraintFunction(obj)
            x = obj.designVariable.value;
            l = obj.lambda;
            N = obj.nElem;
            fx = [x(N+1)-l(1),x(N+1)-l(2),(1/N)*sum(x(1:N))-1]';
            obj.value = fx;
        end
        
      function [V,D] = computeEigenFunctionAndValues(obj,B,K)
            [V,D]=eigs(B,K,2,'SM');
        end        
        
        function [dfdx, dfdx2] = computeConstraintDerivative(obj)
                x = obj.designVariable.value;           
                Belem = obj.Be;
                N = obj.nElem;
                m = obj.nConstraints;
                dfdx=zeros(m,N+1);
                dfdx(3,1:N)=(1/N)*ones(1,N);
                dfdx2 = 0*dfdx;
                if abs(obj.D(2,2)-obj.D(1,1))> 1 % Simple eigenvalues
                    W=zeros(2*N+2,2);
                    for i=3:2*N
                        W(i,1)=obj.v1(i-2);
                    end
                    for i=1:N
                        dfdx(1,i)= -(2*x(i,1))*(W(2*(i-1)+1: 2*(i-1)+4,1)'*Belem*W(2*(i-1)+1: 2*(i-1)+4,1));
                    end
                    for i=3:2*N
                        W(i,2)=obj.v2(i-2);
                    end
                    for i=1:N
                        dfdx(2,i)= -(2*x(i,1))*(W(2*(i-1)+1: 2*(i-1)+4,2)'*Belem*W(2*(i-1)+1: 2*(i-1)+4,2));
                    end  
                else % Dobles eigenvalues
                    obj.D
                    disp('dobles')
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
                    % DERIVATIVES MATRIX DEFINITION
                    for i=1:N
                        dQ1(i,1)= (2*x(i,1))*(Q1(2*(i-1)+1: 2*(i-1)+4,1)'*Belem*Q1(2*(i-1)+1: 2*(i-1)+4,1));
                        dQ2(i,1)= (2*x(i,1))*(Q2(2*(i-1)+1: 2*(i-1)+4,1)'*Belem*Q2(2*(i-1)+1: 2*(i-1)+4,1));
                        dQ1Q2(i,1)= (2*x(i,1))*(Q1(2*(i-1)+1: 2*(i-1)+4,1)'*Belem*Q2(2*(i-1)+1: 2*(i-1)+4,1));
                        A=[dQ1(i,1) dQ1Q2(i,1); dQ1Q2(i,1) dQ2(i,1)];
                        [U,R]=eigs(A,2,'SM');
                        S=sort(diag(R));
                        dfdx(1,i)=-S(1);
                        dfdx(2,i)=-S(2);
                    end
                end  
                dfdx(1,N+1)=1;
                dfdx(2,N+1)=1; 
                obj.gradient = dfdx;
        end  
        
        function [v1,v2] = reorderModes(obj,lambda,V,D)
            if lambda(1)==D(1,1)
                v1=V(:,1);
                v2=V(:,2);
            else
                v1=V(:,2);
                v2=V(:,1);
            end
        end
        
       function computeBucklingModes(obj)
            N = obj.nElem;
            Mode1=zeros(2*N+2);
            Mode2=zeros(2*N+2);
            for i=3:2*N
                Mode1(i)=obj.v1(i-2);
                Mode2(i)=obj.v2(i-2);
            end
            obj.mode1 = Mode1;
            obj.mode2 = Mode2;
        end        
        
        
    end
    
end