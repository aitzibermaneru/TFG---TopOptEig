classdef EigModes < handle
    
    properties (Access = public)
        lambda
        V
        v1
        v2
        D
    end
    
    properties (Access = private)
        nElem
        freeNodes
        length
        youngModulus
        inertiaMoment
        stiffnessMatrix
        bendingMatrix
        elementalBendingMatrix
        e
        E1
        E2
        mode1
        mode2
    end

    properties (Access = private)
        designVariable
        stiffnessMatComputer
        bendingMatComputer
    end

    methods (Access = public)
        
        function obj = EigModes(cParams)
            obj.init(cParams)
            obj.createStiffnessMatrix();
            obj.createBendingMatrix();
        end
        
        function compute(obj,nIter)
            obj.computeEigenModesAndValues();
            obj.computeBucklingModes();
            obj.computeConvergence2Eigenvalues(nIter);
        end

        function plot(obj,x)
            N = obj.nElem;
            L = obj.length;
            % axis and profile
            ch = 0:L:1-L;
            h  = 0:L:1;
            z = sqrt(x(1:N));
            % Plotting Clamped-clamped configuration
            figure(1)
            subplot(2,2,[1 3]);plot(ch,z)
            title('Clamped-Clamped Column Profile','Interpreter', 'latex','FontSize',20, 'fontweight','b');
            xlabel('x','Interpreter', 'latex','fontsize',14,'fontweight','b');
            ylabel('A(x)','Interpreter', 'latex','fontsize',14,'fontweight','b');
            % Buckling modes
            subplot(2,2,2); plot(h,-obj.mode1(1:2:2*N+2));
            title('First Buckling Mode','Interpreter', 'latex','FontSize',14, 'fontweight','b')
            subplot(2,2,4); plot(h,-obj.mode2(1:2:2*N+2));
            title('Second Buckling Mode','Interpreter', 'latex','FontSize',14, 'fontweight','b')
            figure(2)
            hold on
            plot(obj.e,obj.E1);
            plot(obj.e,obj.E2);
            hold off
            xlabel('Number of Iteration','Interpreter', 'latex','fontsize',18,'fontweight','b');
            ylabel('Eigenvalues','Interpreter', 'latex','fontsize',18,'fontweight','b');
            axis([0 65 0 100]);
        end

        function fx = provideFunction(obj,iter,eigNum)
            [l] = obj.computeSettings(iter);
            x = obj.designVariable.value;
            N = obj.nElem;
            fx = x(N+1)-l(eigNum);
        end

        function [l] = computeSettings(obj,iter)
            obj.stiffnessMatComputer.compute();
            obj.stiffnessMatrix = obj.stiffnessMatComputer.stiffnessMatrix;
            obj.bendingMatComputer.compute();
            obj.elementalBendingMatrix = obj.bendingMatComputer.elementalBendingMatrix;
            obj.bendingMatrix = obj.bendingMatComputer.bendingMatrix;  
            obj.compute(iter);
            l = obj.lambda;
        end


        function grad = provideDerivative(obj,iter,eigNum)
            Belem = obj.elementalBendingMatrix;
            x = obj.designVariable.value;
            N = obj.nElem;
            if abs(obj.D(2,2)-obj.D(1,1))> 1
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
            else
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
                A  = zeros(2,2);
                for i=1:N
                    dQ1(i,1)= (2*x(i,1))*(Q1(2*(i-1)+1: 2*(i-1)+4,1)'*Belem*Q1(2*(i-1)+1: 2*(i-1)+4,1));
                    dQ2(i,1)= (2*x(i,1))*(Q2(2*(i-1)+1: 2*(i-1)+4,1)'*Belem*Q2(2*(i-1)+1: 2*(i-1)+4,1));
                    dQ1Q2(i,1)= (2*x(i,1))*(Q1(2*(i-1)+1: 2*(i-1)+4,1)'*Belem*Q2(2*(i-1)+1: 2*(i-1)+4,1));
                    A = [dQ1(i,1) dQ1Q2(i,1); dQ1Q2(i,1) dQ2(i,1)];
                    [U,R] = eigs(A,2,'SM');
                    S = sort(diag(R));
                    dfdx(1,i)=-S(1);
                    dfdx(2,i)=-S(2);
                end
            end
            dfdx(1,N+1)=1;
            dfdx(2,N+1)=1;
            grad = dfdx(eigNum,:);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.freeNodes      = cParams.freeNodes;
            obj.nElem          = cParams.nElem;
            obj.length         = cParams.length;
            obj.youngModulus   = cParams.youngModulus;
            obj.inertiaMoment  = cParams.inertiaMoment;
            obj.designVariable = cParams.designVariable;
            obj.e         = zeros(0);
            obj.E1        = zeros(0);
            obj.E2        = zeros(0);
        end

         function createStiffnessMatrix(obj)
            s.nElem          = obj.nElem;
            s.length         = obj.length;
            s.youngModulus   = obj.youngModulus;
            s.inertiaMoment  = obj.inertiaMoment;
            obj.stiffnessMatComputer = StiffnessMatrixComputer(s);

        end

        function createBendingMatrix(obj)
            s.nElem          = obj.nElem;
            s.length         = obj.length;
            s.youngModulus   = obj.youngModulus;
            s.inertiaMoment  = obj.inertiaMoment;
            s.designVariable = obj.designVariable;
            obj.bendingMatComputer = BendingMatrixComputer(s);
        end

        function createConvengeceParam(obj,nIter)
            obj.e  = zeros(nIter);
            obj.E1 = zeros(nIter);
            obj.E2 = zeros(nIter);
        end
        
        function computeEigenModesAndValues(obj)            
            free   = obj.freeNodes;
            Bfree  = obj.bendingMatrix(free,free);
            Kfree  = obj.stiffnessMatrix(free,free);
            [V,D]  = obj.computeEigenFunctionAndValues(Bfree,Kfree);
            lambda = sort(diag(D));
            [v1,v2] = obj.reorderModes(lambda,V,D);
            obj.v1 = v1;
            obj.v2 = v2;
            obj.V  = V;
            obj.D  = D;
            obj.lambda = lambda;
        end

        function [V,D] = computeEigenFunctionAndValues(obj,B,K)
            [V,D] = eigs(B,K,2,'SM');
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
        
        function [v1,v2] = reorderModes(obj,lambda,V,D)
            if lambda(1)==D(1,1)
                v1=V(:,1);
                v2=V(:,2);
            else
                v1=V(:,2);
                v2=V(:,1);
            end
        end

        function computeConvergence2Eigenvalues(obj,iter)
            obj.e(iter)=iter;
            obj.E1(iter)= obj.D(1,1);
            obj.E2(iter)= obj.D(2,2);
        end

    end
    
end