%% ---------- OPTIMAL DESIGN OF A COLUMN AGAINST BUCKLING ------
%
% Eigenvalue optimization solver, focus on an optimal 
% design of a column against buckling 


function EigOptimizationRunner
s.input = loadInputData();
EigOptimizationSolver = EigOptComputer(s);
EigOptimizationSolver.compute();
end

function s = loadInputData()
s.nel   = 200;
s.alpha = 0.25;
s.beta  = 10;
end