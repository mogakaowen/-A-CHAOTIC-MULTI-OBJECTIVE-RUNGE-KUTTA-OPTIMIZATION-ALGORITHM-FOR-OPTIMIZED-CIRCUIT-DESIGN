clc;
% Define problem-specific variables
nP = 50;    % population size
MaxIt = 200; % maximum number of iterations


ChaosVec=zeros(10,MaxIt);
%Calculate chaos vector
for i=1:10
    ChaosVec(i,:)=chaos(i,MaxIt,1);
end

% Define the objective function
Func_name='F21';

[lb,ub,dim,fobj] = BenchmarkFunctions2(Func_name);

results1=[];
results2=[];

% Run the algorithm 10 times
for i = 1:10
    [Best_Cost,Best_X,Convergence_curve] = CMRUN(nP,MaxIt,lb,ub,dim,fobj,ChaosVec(1,:));
    results1(i) = Convergence_curve(1, MaxIt);
    results2(i) = Convergence_curve(2, MaxIt);
    fprintf('Run %d:\n', i);
    fprintf('  Best Cost: %.4f, %.4f\n', Best_Cost);
end

results1'
results2'
