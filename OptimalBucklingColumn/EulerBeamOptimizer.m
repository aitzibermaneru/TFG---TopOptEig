classdef EulerBeamOptimizer < handle
    
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        nElem
        nConstraints
        length
        youngModulus
        inertiaMoment
    end
    
    methods (Access = public)
        
        function obj = EulerBeamOptimizer()
            obj.init()
            obj.bound()
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.nElem        = 100;
            obj.nConstraints = 3;
            obj.length       = 1/obj.nElem;
            obj.youngModulus = 1;
            obj.inertiaMoment = 1;            
        end
        
        function Be = computeBendingMatrix(obj)
            L = obj.length;
            E = obj.youngModulus;
            I = obj.inertiaMoment;
            Be =(E*I/(L^3))*[12 6*L -12 6*L;  6*L 4*L^2 -6*L 2*L^2; -12 -6*L 12 -6*L; 6*L 2*L^2 -6*L 4*L^2];
        end
        
        function Ke = computeStiffnessMatrix(obj)
            L = obj.length;                     
            Ke = 1/(30*L)*[36 3*L -36 3*L; 3*L 4*L^2 -3*L -L^2; -36 -3*L 36 -3*L; 3*L -L^2 -3*L 4*L^2];                        
        end
        
        function [x] = bound(obj)
            N = obj.nElem;
            % This algorithm maximizes the least eigenvalue of a clamped-clamped
            % column, accounting for the presence of simple or multiple eigenvalues.
            % It also illustrates the best strongest profile of that column against
            % buckling as well as it gives its buckling modes. Its input variable is N,
            % which refers to the number of elements of column with the Finite Element
            % Analysis.
            
            
            % CONSTANTS DEFINITION

            L=1/N ;         %Element's length

            
            % VARIABLES DEFINITION
            n_val=N+1;
            m = obj.nConstraints;
            
            % PUNCTUAL LIMITATIONS OF THE DESIGN
            %(Pre-defined in MMA file)
            alpha=0.25;
            beta=10;
            x=ones(N+1,1);
            
            xmin=alpha*ones(N,1);
            xmin=[xmin; 0];
            xmax=beta*ones(N,1);
            xmax=[xmax; 1000];
            xold1=x;
            xold2=x;
            loop=0;
            
            low = zeros(n_val,1);
            upp = ones(n_val,1);
            a0 = 1;
            a_mma = zeros(m,1);
            d = zeros(m,1);
            c = 1000*ones(m,1);
            
            % AUXILIAR VECTORS DEFINITION FOR EIGENVALUES COMPARATIVE
            e=zeros(loop);
            E1=zeros(loop);
            E2=zeros(loop);
            
            Be = obj.computeBendingMatrix();
            Ke = obj.computeStiffnessMatrix();
            
            % ITERATIVE PROCESS
            change = 1;
            while (change > 0.0005) & (loop < 1000)
                loop = loop + 1;
                
                % BENDING AND STIFFNESS MATRICES DEFINTION
                B=sparse(2*N+2, 2*N+2);
                K=sparse(2*N+2, 2*N+2);
                
                % BENDING AND STIFFNESS MATRICES ASSEMBLY
                for elx=1:N
                    edof=[2*elx-1; 2*elx; 2*(elx+1)-1; 2*(elx+1)];
                    B(edof,edof)=B(edof,edof)+(x(elx)^2)*Be;
                    K(edof,edof)=K(edof,edof)+ Ke;
                end
                
                % BOUNDARY CONDITIONS
                fixnodes = union([1,2], [2*N+1,2*N+2]);
                nodes      = 1:2*N+2;
                freenodes= setdiff(nodes,fixnodes);
                
                [V,D,v1,v2,lambda] = obj.computeEigenModesAndValues(freenodes,K,B);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %__________________MMA____________________
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % UPDATING OF THE SOLUTION
                xval = x;
                
                % OBJECTIVE FUNCTION
                [f0val,df0dx,df0dx2] = obj.computeCost(x);

                
                % CONSTRAINTS VECTOR
                fval=obj.computeConstraintFunction(x,lambda);
                
                [dfdx,dfdx2] = obj.computeConstraintDerivative(Be,D,x,v1,v2,V);
                
                % INVOKING MMA
%                obj.computeNewDesign()
                [xmma,~,~,~,~,~,~,~,~,low,upp] = ...
                    mmasub(m,n_val,loop,xval,xmin,xmax,xold1,xold2, ...
                    f0val,df0dx,df0dx2,fval,dfdx,dfdx2,low,upp,a0,a_mma,c,d);
                
                
                %CALCULATION THE OF BUCKLING MODES
                Mode1=zeros(2*N+2);
                Mode2=zeros(2*N+2);
                
                for i=3:2*N
                    Mode1(i)=v1(i-2);
                    Mode2(i)=v2(i-2);
                end
                
                % COLUMN'S PROFILE
                z=sqrt(x(1:N));

                % CONVERGENCE OF DOUBLE EIGENVALUES
                e(loop)=loop;
                E1(loop)= D(1,1);
                E2(loop)=D(2,2);                
                
                obj.plot(z,Mode1,Mode2,e,E1,E2)
                %OUTPUT VARIABLES UPDATING
                xold2 = xold1;
                xold1 = xval;
                
                x = xmma;
                change = max(abs(x-xold1));
                
                
                % PRINTING OF THE RESULTS
                obj.displayIteration(loop,x,f0val,D)

                cost(loop) = -xmma(N+1);
                vol(loop) = (1/N)*sum(x(1:N));
                figure(3)
                plot(cost)
                figure(4)
                plot(vol)
            end

            
        end

        function [V,D,v1,v2,lambda] = computeEigenModesAndValues(obj,freenodes,K,B)
            [V,D]=eigs(B(freenodes,freenodes),K(freenodes,freenodes),2,'SM');

            lambda=sort(diag(D));

            if lambda(1)==D(1,1)
                v1=V(:,1);
                v2=V(:,2);
            else
                v1=V(:,2);
                v2=V(:,1);
            end

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



        function fx = computeConstraintFunction(obj,x,lambda)
            N = obj.nElem;
            fx = [x(N+1)-lambda(1),x(N+1)-lambda(2),(1/N)*sum(x(1:N))-1]';
        end
        
        function [dfdx, dfdx2] = computeConstraintDerivative(obj,Be,D,x,v1,v2,V)
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
                        dfdx(1,i)= -(2*x(i,1))*(W(2*(i-1)+1: 2*(i-1)+4,1)'*Be*W(2*(i-1)+1: 2*(i-1)+4,1));
                    end
                    
                    for i=3:2*N
                        W(i,2)=v2(i-2);
                    end
                    
                    for i=1:N
                        dfdx(2,i)= -(2*x(i,1))*(W(2*(i-1)+1: 2*(i-1)+4,2)'*Be*W(2*(i-1)+1: 2*(i-1)+4,2));
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
                        dQ1(i,1)= (2*x(i,1))*(Q1(2*(i-1)+1: 2*(i-1)+4,1)'*Be*Q1(2*(i-1)+1: 2*(i-1)+4,1));
                        dQ2(i,1)= (2*x(i,1))*(Q2(2*(i-1)+1: 2*(i-1)+4,1)'*Be*Q2(2*(i-1)+1: 2*(i-1)+4,1));
                        dQ1Q2(i,1)= (2*x(i,1))*(Q1(2*(i-1)+1: 2*(i-1)+4,1)'*Be*Q2(2*(i-1)+1: 2*(i-1)+4,1));
                        
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
        
        function computeNewDesign()
            
        end
        
        function displayIteration(obj,loop,x,f0val,D)
            N = obj.nElem;
            disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',f0val) ...
                ' Vol.: ' sprintf('%6.3f',  (1/N)*(sum(x)-x(N+1))  ) ...
                ' ch.: ' sprintf('%6.3f',abs(D(2,2)-D(1,1)) )])
            
        end
        
         function plot(obj,z,Mode1,Mode2,e,E1,E2)
                N = obj.nElem;
                L = obj.length;
                % AXES DEFINION FOR FIGURES
                ch= 0:L:1-L;
                h= 0:L:1;
                
                
                % PLOT OF THE BEST STRONGEST PROFILE OF THE COLUMN AGAINST BUCKLING
                % Clamped-clamped configuration
                figure(1)
                
                subplot(2,2,[1 3]);plot(ch,z)
                title('Clamped-Clamped Column Profile','Interpreter', 'latex','FontSize',20, 'fontweight','b');
                xlabel('x','Interpreter', 'latex','fontsize',14,'fontweight','b');
                ylabel('A(x)','Interpreter', 'latex','fontsize',14,'fontweight','b');
                
                %Buckling modes
                subplot(2,2,2); plot(h,-Mode1(1:2:2*N+2));
                title('First Buckling Mode','Interpreter', 'latex','FontSize',14, 'fontweight','b')
                
                subplot(2,2,4); plot(h,-Mode2(1:2:2*N+2));
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