classdef IterativeProcessComputer < handle 

     properties (Access = public)
         designVariable
%         dualVariable
%         cost
%         constraint
        % optimizer
         optimizerType
%         incrementalScheme
%         optimizerSettings
     end
    
    properties (Access = private)
        nElem
        nConstraints
        length
        nValues
        youngModulus
        inertiaMoment
        minThick
        maxThick
        mode1
        mode2
        freeNodes
        v1
        v2
        V
        lambda
        D
        bendingMatrix
        stiffnessMatrix
        Be
        nIter
        maxIter
        e
        E1
        E2
    end

    properties (Access = private)
        xMin
        xMax
        xOld1
        xOld2
        lOW
        uPP
        a0Val
        aMMA
        dVal
        cVal
        hasFinished
        change
        Cost
        cost
        constraint
        xMMA
        xVal
        vol
    end

    methods (Access = public)

        function obj = IterativeProcessComputer(cParams)
            obj.init(cParams);
%             obj.createIncrementalScheme(cParams);
%             obj.createDesignVariable(cParams);
%             obj.createHomogenizedVarComputer(cParams)
%             obj.createCostAndConstraint(cParams);
%             obj.createDualVariable();
%             obj.createOptimizer(cParams);
%             obj.createVideoMaker(cParams);
        end

        function compute(obj)
            obj.computeIterativeProcess();
        end
    end

    methods (Access = private)

         function init(obj,cParams)
            obj.nElem          = cParams.nElem;
            obj.nConstraints   = cParams.nConstraints;
            obj.length         = cParams.length;
            obj.youngModulus   = cParams.youngModulus;
            obj.inertiaMoment  = cParams.inertiaMoment;
            obj.minThick       = cParams.minThick;
            obj.maxThick       = cParams.maxThick;
            obj.nValues        = cParams.nValues;
            obj.nIter          = cParams.loop;
            obj.designVariable = cParams.designVariable;
            obj.xMin           = cParams.mmaParams.xMin;
            obj.xMax           = cParams.mmaParams.xMax;
            obj.xOld1          = cParams.mmaParams.xOld1;
            obj.xOld2          = cParams.mmaParams.xOld2;
            obj.lOW            = cParams.mmaParams.lOW;
            obj.uPP            = cParams.mmaParams.uPP;
            obj.a0Val          = cParams.mmaParams.a0Val;
            obj.aMMA           = cParams.mmaParams.aMMA;
            obj.dVal           = cParams.mmaParams.dVal;
            obj.cVal           = cParams.mmaParams.cVal;
            obj.maxIter        = cParams.maxIter;
            obj.optimizerType  = cParams.optimizerType;
          % obj.designVariable = 'BucklingArea';
         end

         function obj = computeIterativeProcess(obj)
             % Step 1
             % create cost and constraint classes
             % Step 2
%             s.type = 'MMA';
%             s.incScheme = obj.createIncremenalScheme();
%             s.designVar = obj.createDesignVariable();
%             s.homogenizedVarComputer = obj.HomogenizedVarComputer();
%             s.costConstraint = obj.createCostConstraint();
%             s.dualVariable = obj.createDualVariable();
%             s.optimizer = obj.createOptimizer();
%             opt = Optimizer(s);
%             obj.optimizer.solveProblem();

            % s.designVar = obj.designVariable;
            % s.type = obj.optimizerType;
%             s.constraintCase = ;
%             s.cost = ;
%             s.constraint = ;
%             s.dualVariable = ;
%             s.maxIter = ;
%             s.incrementalScheme = ;
%             s.targetParameters = ;
%             s.historyPrinterSettings = ;
            % obj.optimizer = Optimizer.create();

             obj.updateSetting();
             obj.change = 1;
             obj.hasFinished = 0;
             while ~obj.hasFinished
                obj.increaseIter();
                obj.update();
                obj.updateStatus();
                obj.computeNewDesign();
                obj.updateOutput();
                obj.displayIteration()
                obj.plotFigures();
            end

        end

        function updateSetting(obj)
            obj.e  = zeros(obj.nIter);
            obj.E1 = zeros(obj.nIter);
            obj.E2 = zeros(obj.nIter);
            obj.computeElementalBendingMatrix();
            obj.computeStiffnessMatrix();
            obj.computeBoundaryConditions();
        end

        function increaseIter(obj)
            obj.nIter = obj.nIter+1;
        end

        function update(obj)
            obj.computeBendingMatrix();
            obj.computeEigenModesAndValues();
            obj.computeBucklingModes();
            obj.computeConvergence2Eigenvalues();
        end

        function updateStatus(obj)
            obj.hasFinished = (obj.change <= 0.0005) || (obj.nIter >= obj.maxIter);
        end

        function updateOutput(obj)
            obj.xOld2 = obj.xOld1;
            obj.xOld1 = obj.xVal;
            obj.designVariable = obj.xMMA;
            obj.change = max(abs(obj.designVariable-obj.xOld1));
        end

        function plotFigures(obj)
            N = obj.nElem;
            iter = obj.nIter;
            x = obj.designVariable;
            obj.cost.val(iter) = -obj.xMMA(N+1);
            obj.vol(iter) = (1/N)*sum(x(1:N));
            obj.plot()
            figure(3)
            plot(obj.cost.val)
            figure(4)
            plot(obj.vol)
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

        function computeBendingMatrix(obj)
            x = obj.designVariable;
            N = obj.nElem;
            B=sparse(2*N+2, 2*N+2);    
            for iElem = 1: N
                iDof=[2*iElem-1; 2*iElem; 2*(iElem+1)-1; 2*(iElem+1)];
                B(iDof,iDof)=B(iDof,iDof)+(x(iElem)^2)*obj.Be;
            end
            obj.bendingMatrix = B;
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

        function  computeBoundaryConditions(obj)
            N = obj.nElem;
            fixnodes = union([1,2], [2*N+1,2*N+2]);
            nodes = 1:2*N+2;
            free  = setdiff(nodes,fixnodes);
            obj.freeNodes = free;
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

        function [v1,v2] = reorderModes(obj,lambda,V,D)
            if lambda(1)==D(1,1)
                v1=V(:,1);
                v2=V(:,2);
            else
                v1=V(:,2);
                v2=V(:,1);
            end
        end

        function [V,D] = computeEigenFunctionAndValues(obj,B,K)
            [V,D]=eigs(B,K,2,'SM');
        end

        function computeNewDesign(obj)
            n_val=obj.nValues;
            m = obj.nConstraints;
            xmin = obj.xMin;
            xmax = obj.xMax;
            xold1 = obj.xOld1;
            xold2 = obj.xOld2;
            low = obj.lOW;
            upp = obj.uPP;
            a0 = obj.a0Val;
            a_mma = obj.aMMA;
            d = obj.dVal;
            c = obj.cVal;
            iter = obj.nIter;
            x = obj.designVariable;
            xval = x; % design variable -----> será el input al optimizier_MMA

            obj.computeCost(); 
            fval =obj.computeConstraintFunction(x);
            [dfdx,dfdx2] = obj.computeConstraintDerivative(x);

            f0val = obj.cost.value;
            df0dx = obj.cost.gradient;
            [xmma,~,~,~,~,~,~,~,~,low,upp] = ...
                mmasub(m,n_val,iter,xval,xmin,xmax,xold1,xold2, ...
                f0val,df0dx,df0dx2,fval,dfdx,dfdx2,low,upp,a0,a_mma,c,d);
            obj.lOW = low;
            obj.uPP = upp;
            obj.xOld1 = xold1;
            obj.xVal = xval;
            obj.xMMA = xmma;
        end

        function computeCost(obj)
%             N = obj.nElem;
%             f0val = -x(N+1); 
%             df0dx = zeros(N+1,1);
%             df0dx(N+1) = -1;
%             df0dx2 = 0*df0dx;
%             obj.cost.value = f0val;
%             obj.cost.gradient = df0dx;
              s.nElem = obj.nElem;
              s.designVariable = obj.designVariable;
              solution = Cost(s);
        end

        function fx = computeConstraintFunction(obj,x)
            l = obj.lambda;
            N = obj.nElem;
            fx = [x(N+1)-l(1),x(N+1)-l(2),(1/N)*sum(x(1:N))-1]';
            obj.constraint.value = fx;
        end
        
        function [dfdx, dfdx2] = computeConstraintDerivative(obj,x)
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
                obj.constraint.gradient = dfdx;
        end
        
        function displayIteration(obj)
            x = obj.designVariable;
            N = obj.nElem;
            f0val = obj.cost.value;
            iter = obj.nIter;
            disp([' It.: ' sprintf('%4i',iter) ' Obj.: ' sprintf('%10.4f',f0val) ...
                ' Vol.: ' sprintf('%6.3f',  (1/N)*(sum(x)-x(N+1))  ) ...
                ' ch.: ' sprintf('%6.3f',abs(obj.D(2,2)-obj.D(1,1)) )])
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

        function computeConvergence2Eigenvalues(obj)
            iter = obj.nIter;
            obj.e(iter)=iter;
            obj.E1(iter)= obj.D(1,1);
            obj.E2(iter)= obj.D(2,2);
        end

        function plot(obj)
                x = obj.designVariable;
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

end