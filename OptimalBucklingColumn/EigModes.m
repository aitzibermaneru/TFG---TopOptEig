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
        stiffnessMatrix
        bendingMatrix
        e
        E1
        E2
        mode1
        mode2
    end
    
    methods (Access = public)
        
        function obj = EigModes(cParams)
            obj.init(cParams)
        end
        
        function compute(obj,K,B,nIter)
            obj.stiffnessMatrix = K;
            obj.bendingMatrix = B;
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
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.freeNodes = cParams.freeNodes;
            obj.nElem     = cParams.nElem;
            obj.length    = cParams.length;
            obj.e         = zeros(0);
            obj.E1        = zeros(0);
            obj.E2        = zeros(0);
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