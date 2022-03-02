classdef Constraint < handle
    
    properties (GetAccess = public, SetAccess = private)
        mode1
        mode2
        D
        gradient
        value
    end

    properties (Access = private) % classes
        designVariable
        bendingMat
        stiffnessMat
        eigModes
        constraintDeriv
    end
    
    properties (Access = private) % computed
        bendingMatrix
        elementalBendingMatrix
        stiffnessMatrix         
        v1
        v2
        V
        lambda
    end
    
    properties (Access = private) % input  
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
            obj.createStiffnessMatrix();   
            obj.createBendingMatrix();
            obj.createEigModes();
            obj.createConstraintDerivative();
        end
        
        function computeFunctionAndGradient(obj)
            obj.computeBendingMatrix();
            obj.computeEigModes();           
            obj.computeConstraintFunction();
            obj.computeConstraintDerivative();
            %obj.eigenValueModes.plot();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.designVariable = cParams.designVariable;
            obj.nElem          = cParams.nElem;
            obj.freeNodes      = cParams.freeNodes;
            obj.length         = cParams.length;
            obj.youngModulus   = cParams.youngModulus;
            obj.inertiaMoment  = cParams.inertiaMoment;   
            obj.nConstraints   = cParams.nConstraints;
        end

        function createStiffnessMatrix(obj)
            s.nElem          = obj.nElem;
            s.length         = obj.length;
            s.youngModulus   = obj.youngModulus;
            s.inertiaMoment  = obj.inertiaMoment;
            obj.stiffnessMat = StiffnessMatrixComputer(s);
            obj.stiffnessMat.compute();
            obj.stiffnessMatrix = obj.stiffnessMat.stiffnessMatrix;
        end

        function createBendingMatrix(obj)
            s.nElem          = obj.nElem;
            s.length         = obj.length;
            s.youngModulus   = obj.youngModulus;
            s.inertiaMoment  = obj.inertiaMoment;
            s.designVariable = obj.designVariable;
            obj.bendingMat = BendingMatrixComputer(s);
        end 

        function createEigModes(obj)
            s.freeNodes  = obj.freeNodes;
            s.nElem      = obj.nElem;
            obj.eigModes = EigModes(s);
        end

        function createConstraintDerivative(obj)
            s.nElem = obj.nElem;
            s.m     = obj.nConstraints;
            obj.constraintDeriv = ConstraintDerivative(s);
        end


        function computeBendingMatrix(obj)
           obj.bendingMat.compute();
           obj.bendingMatrix =  obj.bendingMat.bendingMatrix;
           obj.elementalBendingMatrix = obj.bendingMat.elementalBendingMatrix;
        end 

        function computeEigModes(obj)
            K = obj.stiffnessMatrix;
            B = obj.bendingMatrix;
            obj.eigModes.compute(K,B);
            obj.lambda = obj.eigModes.lambda; 
            obj.D      = obj.eigModes.D;
            obj.V      = obj.eigModes.V;
            obj.v1     = obj.eigModes.v1;
            obj.v2     = obj.eigModes.v2;
            obj.mode1  = obj.eigModes.mode1;
            obj.mode2  = obj.eigModes.mode2;
        end

        function fx = computeConstraintFunction(obj)
            x = obj.designVariable.value;
            l = obj.lambda;
            N = obj.nElem;
            fx = [x(N+1)-l(1),x(N+1)-l(2),(1/N)*sum(x(1:N))-1]';
            obj.value = fx;
        end

        function [dfdx, dfdx2] = computeConstraintDerivative(obj)
            x = obj.designVariable.value;
            Belem = obj.elementalBendingMatrix;
            if abs(obj.D(2,2)-obj.D(1,1))> 1 % Simple eigenvalues
                obj.constraintDeriv.computeSimple(x,Belem,obj.v1,obj.v2);
            else % Dobles eigenvalues
                obj.constraintDeriv.computeDouble(x,Belem,obj.V);
            end
            dfdx  = obj.constraintDeriv.dfdx;
            dfdx2 = obj.constraintDeriv.dfdx2;
            dfdx(1,obj.nElem+1)=1;
            dfdx(2,obj.nElem+1)=1;
            obj.gradient = dfdx';
        end

    end
    
end