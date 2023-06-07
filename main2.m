clear 
close all
clc

nP = 100;                           % Number of population
MaxIt = 1000;                       % Maximum number of iterations
lb = 0;                           % Lower bound of decision variables
ub = 100;                        % Upper bound of decision variables
dim = 50;                           % Number of decision variables

ChaosVec=zeros(10,MaxIt);
%Calculate chaos vector
for i=1:10
    ChaosVec(i,:)=chaos(i,MaxIt,1);
end

fobj=@circuit;

[Best_Cost,Best_X,Convergence_curve] = CMRUN(nP,MaxIt,lb,ub,dim,fobj,ChaosVec(1,:));

%% Draw objective space

%Pareto Front
figure()
plot(Convergence_curve(1,:), Convergence_curve(2,:), 'ko')
title('Pareto Front')
xlabel('Noise');
ylabel('Distortion');


weights = [0.5,0.5];
Fitness = ones(1,MaxIt);
for i=1:MaxIt
    Fitness(1,:) = weights(1,1)*Convergence_curve(1,:) + weights(1,2)*Convergence_curve(2,:);
end



%Pareto Sum Fitness
figure()
plot(Fitness,'Color','b','LineWidth',2)
title('Chaotic Multi-Objective RUN Optimization Algorithm')
xlabel('Iterations');
ylabel('Pareto Sum Fitness');
grid on
box on
legend('CMRUN')
