%RUN and CRUN main Code
%Owen Mogaka Nyandieka, Department of Electrical and Information Engineering University Of Nairobi
clear 
close all
clc

nP=50;          % Number of Population

Func_name='F1'; % Name of the test function, range from F1-F14

MaxIt=500;      % Maximum number of iterations

% Load details of the selected benchmark function
[lb,ub,dim,fobj]=BenchmarkFunctions(Func_name);

[Best_fitness,BestPositions,Convergence_curve] = CRUN(nP,MaxIt,lb,ub,dim,fobj);

%% Draw objective space

figure,
hold on
semilogy(Convergence_curve,'Color','r','LineWidth',4);
title('Convergence curve')
xlabel('Iteration');
ylabel('Best fitness obtained so far');
axis tight
grid off
box on
legend('RUN')


