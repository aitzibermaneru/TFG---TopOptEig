classdef Sh_doubleFirstEig < ShapeFunctional
    
    properties (Access = private)
        nElem
        designVariable
    end
    
    properties (Access = private)
        
    end   
    
    methods (Access = public)
        
        function obj = Sh_doubleFirstEig(cParams)
            obj.init(cParams)
        end

        function computeFunctionAndGradient(obj)
            obj.computeFunction();
            obj.computeGradient();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.nElem = cParams.nElem;
            obj.designVariable = cParams.designVariable;
        end

        function computeFunction(obj)
            fx = x(N+1)-l(1);
            obj.value = fx;
        end

        function computeGradient(obj)
            if abs(obj.D(2,2)-obj.D(1,1))> 1
                for i=3:2*N
                    W(i,1)=obj.v1(i-2);
                end
                for i=1:N
                    dfdx(1,i)= -(2*x(i,1))*(W(2*(i-1)+1: 2*(i-1)+4,1)'*Belem*W(2*(i-1)+1: 2*(i-1)+4,1));
                end
            else
                disp('dobles')
                Q1=zeros(2*N+2,1);
                Q2=zeros(2*N+2,1);
                dQ1=zeros(N,1);
                dQ2=zeros(N,1);
                dQ1Q2=zeros(N,1);
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