classdef IterativeProcessComputer < handle 

    properties (Access = public)
        designVariable
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
    end

    properties (Access = private)
        xMin
        xMax
        xOld1
        xOld2
        looping
        lOW
        uPP
        a0Val
        aMMA
        dVal
        cVal
    end

    methods (Access = public)

        function obj = IterativeProcessComputer(cParams)
            obj.init(cParams);
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
            obj.looping        = cParams.loop;
            obj.designVariable = cParams.designVariable;
            obj.xMin           = cParams.xMin;
            obj.xMax           = cParams.xMax;
            obj.xOld1          = cParams.xOld1;
            obj.xOld2          = cParams.xOld2;
            obj.lOW            = cParams.lOW;
            obj.uPP            = cParams.uPP;
            obj.a0Val          = cParams.a0Val;
            obj.aMMA           = cParams.aMMA;
            obj.dVal           = cParams.dVal;
            obj.cVal           = cParams.cVal;
         end

        function obj = computeIterativeProcess(obj)
            N = obj.nElem;
            n_val=obj.nValues;
            m = obj.nConstraints;
            iter = obj.looping;
            x = obj.designVariable;
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
            e          = zeros(iter);
            E1             = zeros(iter);
            E2             = zeros(iter);
            % Elemental Matrices
            obj.computeElementalBendingMatrix();
            obj.computeStiffnessMatrix();

            obj.computeBoundaryConditions();

            % looping
            change =1;
            while (change > 0.0005) && (iter < 1000)
                iter = iter + 1;
                obj.computeBendingMatrix(x);
                obj.computeEigenModesAndValues();
                % MMA
                [xmma,low,upp,xold1,xval,f0val] = obj.computeNewDesign(x,m,n_val,iter,xmin,xmax,xold1,xold2,low,upp,a0,a_mma,c,d);
                % BUCKLING MODES
                obj.computeBucklingModes();

                % CONVERGENCE DOUBLES EIGENVLAUES
                [e,E1,E2]= obj.computeConvergence2Eigenvalues(iter,e,E1,E2);
                % PLOTTING
                %OUTPUT VARIABLES UPDATING
                xold2 = xold1;
                xold1 = xval;
                x = xmma;
                change = max(abs(x-xold1));
                % PRINTING OF THE RESULTS
                obj.displayIteration(iter,x,f0val)
                % PLOTTING, COST AND VOLUME
                cost(iter) = -xmma(N+1);
                vol(iter) = (1/N)*sum(x(1:N));  % x=xmma, esta bien?
                
                obj.plotFigures(e,E1,E2,x,cost,vol)
            end
        end

        

        function plotFigures(obj,e,E1,E2,x,cost,vol)
                obj.plot(e,E1,E2,x)                
                figure(3)
                plot(cost)
                figure(4)
                plot(vol)
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

        function computeBendingMatrix(obj,x)
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

        function [xmma,low,upp,xold1,xval,f0val] = computeNewDesign(obj,x,m,n_val,loop,xmin,xmax,xold1,xold2,low,upp,a0,a_mma,c,d)
%             s.designVar = x;
%             s.type = obj.optimizerType;
%             s.constraintCase = ;
%             s.cost = ;
%             s.constraint = ;
%             s.dualVariable = ;
%             s.maxIter = ;
%             s.incrementalScheme = ;
%             s.targetParameters = ;
%             s.historyPrinterSettings = ;
%             create = Optimizer.create(s);
%             solution = create.solveProblem();
            % UPDATING OF THE SOLUTION
            xval = x; % design variable -----> será el input al optimizier_MMA
            % OBJECTIVE FUNCTION
            [f0val,df0dx,df0dx2] = obj.computeCost(x);  % ya lo calcula lo otro
            % CONSTRAINTS
            fval=obj.computeConstraintFunction(x);  % ya lo calcula
            [dfdx,dfdx2] = obj.computeConstraintDerivative(x); % ya lo calcula
            % INVOKING MMA
            [xmma,~,~,~,~,~,~,~,~,low,upp] = ...
                mmasub(m,n_val,loop,xval,xmin,xmax,xold1,xold2, ...
                f0val,df0dx,df0dx2,fval,dfdx,dfdx2,low,upp,a0,a_mma,c,d);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end

        function [f0val,df0dx,df0dx2] = computeCost(obj,x)
            N = obj.nElem;
            f0val=-x(N+1);
            % OBJECTIVE FUNCTION'S FIRST DERIVATIVE
            df0dx=zeros(N+1,1);
            df0dx(N+1)=-1;
            % OBJECTIVE FUNCTION'S SECOND DERIVATIVE
            df0dx2 = 0*df0dx;
        end

        function fx = computeConstraintFunction(obj,x)
            l = obj.lambda;
            N = obj.nElem;
            fx = [x(N+1)-l(1),x(N+1)-l(2),(1/N)*sum(x(1:N))-1]';
        end
        
        function [dfdx, dfdx2] = computeConstraintDerivative(obj,x)
                D = obj.D;
                V = obj.V;
                v1 = obj.v1;
                v2 = obj.v2;
                Belem = obj.Be;
                
                N = obj.nElem;
                m = obj.nConstraints;
                % CONSTRAINTS VECTOR'S FIRST DERIVATIVE
                dfdx=zeros(m,N+1);
                dfdx(3,1:N)=(1/N)*ones(1,N);
                % CONSTRAINTS VECTOR'S SECOND DERIVATIVE
                dfdx2 = 0*dfdx;
                % MULTIPLE EIGENVALUE'S IDENTIFICATION
                if abs(D(2,2)-D(1,1))> 1
                    % OBJECTIVE FUNCTION'S FIRST DERIVATIVE CALCULATION FOR SIMPLE EIGENVALUES
                    W=zeros(2*N+2,2);
                    for i=3:2*N
                        W(i,1)=v1(i-2);
                    end
                    for i=1:N
                        dfdx(1,i)= -(2*x(i,1))*(W(2*(i-1)+1: 2*(i-1)+4,1)'*Belem*W(2*(i-1)+1: 2*(i-1)+4,1));
                    end
                    for i=3:2*N
                        W(i,2)=v2(i-2);
                    end
                    for i=1:N
                        dfdx(2,i)= -(2*x(i,1))*(W(2*(i-1)+1: 2*(i-1)+4,2)'*Belem*W(2*(i-1)+1: 2*(i-1)+4,2));
                    end  
                else
                    D
                    disp('dobles')
                    % OBJECTIVE FUNCTION'S FIRST DERIVATIVE CALCULATION FOR DOUBLE EIGENVALUES
                    % AUXILIAR VECTORS FOR DERIVATIVE'S CALCULATION
                    Q1=zeros(2*N+2,1);
                    Q2=zeros(2*N+2,1);
                    dQ1=zeros(N,1);
                    dQ2=zeros(N,1);
                    dQ1Q2=zeros(N,1);
                    for i=3:2*N
                        Q1(i,1)=V(i-2,1);
                    end
                    for i=3:2*N
                        Q2(i,1)=V(i-2,2);
                    end
                    % DERIVATIVES MATRIX DEFINITION
                    A=zeros(2,2); 
                    for i=1:N
                        %Derivadas.
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
        end
        
        function displayIteration(obj,loop,x,f0val)
            N = obj.nElem;
            disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',f0val) ...
                ' Vol.: ' sprintf('%6.3f',  (1/N)*(sum(x)-x(N+1))  ) ...
                ' ch.: ' sprintf('%6.3f',abs(obj.D(2,2)-obj.D(1,1)) )])
        end

        function computeBucklingModes(obj)
            N = obj.nElem;
            %CALCULATION THE OF BUCKLING MODES
            Mode1=zeros(2*N+2);
            Mode2=zeros(2*N+2);
            for i=3:2*N
                Mode1(i)=obj.v1(i-2);
                Mode2(i)=obj.v2(i-2);
            end
            obj.mode1 = Mode1;
            obj.mode2 = Mode2;
        end

        function [e,E1,E2]= computeConvergence2Eigenvalues(obj,iter,e,E1,E2)
            % CONVERGENCE OF DOUBLE EIGENVALUES
            e(iter)=iter;
            E1(iter)= obj.D(1,1);
            E2(iter)= obj.D(2,2);
        end

        function plot(obj,e,E1,E2,x)
                N = obj.nElem;
                L = obj.length;
                % AXES DEFINION FOR FIGURES
                ch= 0:L:1-L;
                h= 0:L:1;       
                % COLUMN'S PROFILE
                z=sqrt(x(1:N));
                % PLOT OF THE BEST STRONGEST PROFILE OF THE COLUMN AGAINST BUCKLING
                % Clamped-clamped configuration
                figure(1)
                subplot(2,2,[1 3]);plot(ch,z)
                title('Clamped-Clamped Column Profile','Interpreter', 'latex','FontSize',20, 'fontweight','b');
                xlabel('x','Interpreter', 'latex','fontsize',14,'fontweight','b');
                ylabel('A(x)','Interpreter', 'latex','fontsize',14,'fontweight','b');
                %Buckling modes
                subplot(2,2,2); plot(h,-obj.mode1(1:2:2*N+2));
                title('First Buckling Mode','Interpreter', 'latex','FontSize',14, 'fontweight','b')
                subplot(2,2,4); plot(h,-obj.mode2(1:2:2*N+2));
                title('Second Buckling Mode','Interpreter', 'latex','FontSize',14, 'fontweight','b')
                figure(2)
                hold on                
                plot(e,E1);
                plot(e,E2);
                hold off
                xlabel('Number of Iteration','Interpreter', 'latex','fontsize',18,'fontweight','b');
                ylabel('Eigenvalues','Interpreter', 'latex','fontsize',18,'fontweight','b');
                axis([0 65 0 100]);                             
            end        
    end

end