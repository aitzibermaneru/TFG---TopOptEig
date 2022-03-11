classdef Sh_doubleSecondEig < ShapeFunctional
    
    properties (Access = private)
        nElem
        designVariable
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
    end
    
    methods (Access = public)

        function computeFunction(obj)
            fx = x(N+1)-l(2);
            obj.value = fx;
        end

        function computeGradient(obj)

            if abs(obj.D(2,2)-obj.D(1,1))> 1 % Simple eigenvalues
                W=zeros(2*N+2,2);
                for i=3:2*N
                    W(i,2)=obj.v2(i-2);
                end
                for i=1:N
                    dfdx(2,i)= -(2*x(i,1))*(W(2*(i-1)+1: 2*(i-1)+4,2)'*Belem*W(2*(i-1)+1: 2*(i-1)+4,2));
                end
            else % Dobles eigenvalues
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
                        dfdx(2,i)=-S(2);
                    end
                    dfdx(2,N+1)=1;
                    obj.gradient = dfdx';
            end
        end
    end
    
end