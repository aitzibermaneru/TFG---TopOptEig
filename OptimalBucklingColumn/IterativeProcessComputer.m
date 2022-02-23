classdef IterativeProcessComputer < handle 

     properties (Access = public)
         designVariable
%         dualVariable
%         cost
%         constraint
%         optimizer
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
        nIter
        maxIter
        e
        E1
        E2
        freeNodes
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
        costHistory
        cost
        constraint
        xMMA
        xVal
        vol
    end

    methods (Access = public)

        function obj = IterativeProcessComputer(cParams)
            obj.init(cParams);
            obj.computeBoundaryConditions();            
            obj.createCost();
            obj.createConstraint();

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
% Refactor Constraint
% Refactor mmaParams in Iteriative
% Construct Optimizer;

            % s.designVar = obj.designVariable;
            % s.type     = obj.optimizerType;
%             s.constraintCase = 'INEQUALIY';
%             s.cost = obj.cost;
%             s.constraint = obj.constraint;
%             s.dualVariable = dualVariable;
%             s.maxIter = obj.maxIter;
%             s.incrementalScheme = [];
%             s.targetParameters = [];
%             s.historyPrinterSettings = [];
            % obj.optimizer = Optimizer.create();
             %obj.optimizer.solveProblem();

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
        end

        function increaseIter(obj)
            obj.nIter = obj.nIter+1;
        end

        function update(obj)
            obj.computeConvergence2Eigenvalues();
        end

        function updateStatus(obj)
            obj.hasFinished = (obj.change <= 0.0005) || (obj.nIter >= obj.maxIter);
        end

        function updateOutput(obj)
            obj.xOld2 = obj.xOld1;
            obj.xOld1 = obj.xVal;
            obj.designVariable.update(obj.xMMA);
            obj.change = max(abs(obj.designVariable.value-obj.xOld1));
        end

        function plotFigures(obj)
            N = obj.nElem;
            iter = obj.nIter;
            x = obj.designVariable.value;
            obj.costHistory(iter) = obj.cost.value;
            obj.vol(iter) = (1/N)*sum(x(1:N));
            obj.plot()
            figure(3)
            plot(obj.costHistory)
            figure(4)
            plot(obj.vol)
        end

       function  computeBoundaryConditions(obj)
            N = obj.nElem;
            fixnodes = union([1,2], [2*N+1,2*N+2]);
            nodes = 1:2*N+2;
            free  = setdiff(nodes,fixnodes);
            obj.freeNodes = free;
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
            x = obj.designVariable.value;
            xval = x;

            obj.constraint.computeFunctionAndGradient(); 
            fval = obj.constraint.value;
            dfdx = obj.constraint.gradient;
            dfdx2 = 0;            
            
            obj.cost.computeFunctionAndGradient(); 
            f0val = obj.cost.value;
            df0dx = obj.cost.gradient;
            df0dx2 = 0;

            [xmma,~,~,~,~,~,~,~,~,low,upp] = ...
                mmasub(m,n_val,iter,xval,xmin,xmax,xold1,xold2, ...
                f0val,df0dx,df0dx2,fval,dfdx,dfdx2,low,upp,a0,a_mma,c,d);
            obj.lOW = low;
            obj.uPP = upp;
            obj.xOld1 = xold1;
            obj.xVal = xval;
            obj.xMMA = xmma;
        end

        function createCost(obj)
            s.nElem = obj.nElem;
            s.designVariable = obj.designVariable;
            obj.cost = Cost(s);
            %obj.cost.computeFunctionAndGradient();
        end

        function createConstraint(obj)
            s.designVariable = obj.designVariable;
            s.nElem = obj.nElem;
            s.freeNodes = obj.freeNodes;
            s.length = obj.length;
            s.youngModulus = obj.youngModulus;
            s.inertiaMoment = obj.inertiaMoment; 
            s.nConstraints = obj.nConstraints;
            obj.constraint = Constraint(s);
            obj.constraint.computeFunctionAndGradient();
        end

        function displayIteration(obj)
            x = obj.designVariable.value;
            N = obj.nElem;
            f0val = obj.cost.value;
            iter = obj.nIter;
            disp([' It.: ' sprintf('%4i',iter) ' Obj.: ' sprintf('%10.4f',f0val) ...
                ' Vol.: ' sprintf('%6.3f',  (1/N)*(sum(x)-x(N+1))  ) ...
                ' ch.: ' sprintf('%6.3f',abs(obj.constraint.D(2,2)-obj.constraint.D(1,1)) )])
        end

        function computeConvergence2Eigenvalues(obj)
            iter = obj.nIter;
            obj.e(iter)=iter;
            obj.E1(iter)= obj.constraint.D(1,1);
            obj.E2(iter)= obj.constraint.D(2,2);
        end

        function plot(obj)
                x = obj.designVariable.value;
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
                subplot(2,2,2); plot(h,-obj.constraint.mode1(1:2:2*N+2));
                title('First Buckling Mode','Interpreter', 'latex','FontSize',14, 'fontweight','b')
                subplot(2,2,4); plot(h,-obj.constraint.mode2(1:2:2*N+2));
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