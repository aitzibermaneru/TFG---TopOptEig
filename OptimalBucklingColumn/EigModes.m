classdef EigModes < handle
    
    properties (Access = public)
        lambda
        mode1
        mode2
        V
        D
        v1
        v2
    end
    
    properties (Access = private)
        nElem
        freeNodes
        stiffnessMatrix
        bendingMatrix
    end
    
    methods (Access = public)
        
        function obj = EigModes(cParams)
            obj.init(cParams)
        end
        
        function compute(obj,K,B)
            obj.stiffnessMatrix = K;
            obj.bendingMatrix = B;
            obj.computeEigenModesAndValues();
            obj.computeBucklingModes();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.freeNodes = cParams.freeNodes;
            obj.nElem     = cParams.nElem;
        end
        
        function computeEigenModesAndValues(obj)            
            free  = obj.freeNodes;
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
        function [V,D] = computeEigenFunctionAndValues(obj,B,K)
            [V,D]=eigs(B,K,2,'SM');
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
        
    end
    
end