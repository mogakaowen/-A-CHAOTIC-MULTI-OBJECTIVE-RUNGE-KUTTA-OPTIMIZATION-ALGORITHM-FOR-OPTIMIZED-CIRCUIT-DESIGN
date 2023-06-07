%CMRUN MATLAB Code
%Owen Mogaka Nyandieka, Department of Electrical and Information Engineering University Of Nairobi
function [Best_Cost,Best_X,Convergence_curve]=CMRUN(nP,MaxIt,lb,ub,dim,fobj,ChaosVec)

% Define the number of objectives
nObj = 2;


X=initialization(nP,dim,ub,lb);  % Initialize the set of random solutions 
Cost=zeros(nP,nObj);             % Record the Fitness of all Solutions
Convergence_curve=zeros(nObj,MaxIt);

for i=1:nP
    Cost(i,:) = fobj(X(i,:)); % Calculate the Value of Objective Function
end

cost=sort(Cost);
[~,ind] = min(cost(:,1)); % Determine the Best Solution
[Best_Cost,~] = min(cost); % Determine the Best Cost
Best_X = X(ind,:);

Convergence_curve(1) = Best_Cost(1); % Update the convergence curve with the minimum cost of the first objective
Convergence_curve(2) = Best_Cost(2); % Update the convergence curve with the minimum cost of the second objective

%% Main Loop of RUN 
it=1;%Number of iterations
while it<MaxIt
    it=it+1;
    f=20.*exp(-(12.*(it/MaxIt))); % (Eq.17.6) 
    Xavg = mean(X);               % Determine the Average of Solutions
    SF=2.*(0.5-rand(1,nP)).*f;    % Determine the Adaptive Factor (Eq.17.5)
   
    %Update chaos parameter
    chaotic_param=ChaosVec(it);
    
    for i=1:nP         
            [~,ind_l] = min(cost(:,1));
            lBest = X(ind_l,:);  
            
            [A,B,C]=RndX(nP,i);   % Determine Three Random Indices of Solutions
            [~,ind1] = min(cost([A B C]));
            
            Best_Cost(1) = min(cost(:,1));
            Best_Cost(2) = min(cost(:,2));

            % Determine Delta X (Eqs. 11.1 to 11.3)
            gama = chaotic_param.*(X(i,:)-chaotic_param.*(ub-lb)).*exp(-4*it/MaxIt); % Modified using chaotic map 
            Stp= chaotic_param.*((Best_X-chaotic_param.*Xavg)+gama); % Modified using chaotic map
            DelX = chaotic_param.*(2*chaotic_param-1).*abs(Stp);  % Modified using chaotic map
            
            % Determine Xb and Xw for using in Runge Kutta method
            if cost(i,1)<cost(ind1,1)                
                Xb = X(i,:);
                Xw = X(ind1,:);
            else
                Xb = X(ind1,:);
                Xw = X(i,:);
            end
            SM = RungeKutta(Xb,Xw,DelX);   % Search Mechanism (SM) of RUN based on Runge Kutta Method
                        
            L=rand<0.5;  
            Xc = L.*X(i,:)+(1-L).*X(A,:);  % (Eq. 17.3)
            Xm = L.*Best_X+(1-L).*lBest;   % (Eq. 17.4)           
                        
            vec=[1,-1];
            flag = floor(2*rand(1,dim)+1);
            r=vec(flag);                   % An Interger number 
            
            g = 2*rand;
            mu = 0.5+.1*randn(1,dim);
              
            % Improve the solutions 
            Xnew1 = X(i,:)+r.*Xc.*abs(Xm-X(i,:));
            Xnew2 = Xnew1+rand(1,dim).*SF(i).*(Best_X-X(i,:));
            
            % Update the Solution
            Xnew2 = max(Xnew2,lb);
            Xnew2 = min(Xnew2,ub);
            X(i,:) = Xnew2;
            
            
            % Determine New Solution Based on Runge Kutta Method (Eq.18) 
            if chaotic_param<0.5
                Xnew = (Xc+r.*SF(i).*g.*Xc) + SF(i).*(SM) + mu.*(Xm-Xc);
            else
                Xnew = (Xm+r.*SF(i).*g.*Xm) + SF(i).*(SM)+ mu.*(X(A,:)-X(B,:));
            end  
            
        % Check if solutions go outside the search space and bring them back
        FU=Xnew>ub;FL=Xnew<lb;Xnew=(Xnew.*(~(FU+FL)))+ub.*FU+lb.*FL; 
        CostNew=fobj(Xnew);

        
        if CostNew<cost(i,1)
            X(i,:)=Xnew;
            cost(i,:)=CostNew;
        end
        
        
%% Enhanced solution quality (ESQ)  (Eq. 19)      
        if chaotic_param<0.5
            EXP=exp(-5*chaotic_param*it/MaxIt);
            r = floor(Unifrnd(-1,2,1,1));

            u=2*rand(1,dim); 
            w=Unifrnd(0,2,1,dim).*EXP;              
            
            [A,B,C]=RndX(nP,i);
            Xavg=(X(A,:)+X(B,:)+X(C,:))/3;           %(Eq.19-1)  
            
            beta=rand(1,dim);
            Xnew1 = beta.*(Best_X)+(1-beta).*(Xavg); %(Eq.19-2)
            
            for j=1:dim
                if w(j)<1 
                    Xnew2(j) = Xnew1(j)+r*w(j)*abs((Xnew1(j)-Xavg(j))+randn);
                else
                    Xnew2(j) = (Xnew1(j)-Xavg(j))+r*w(j)*abs((u(j).*Xnew1(j)-Xavg(j))+randn);
                end
            end
            
            FU=Xnew2>ub;FL=Xnew2<lb;Xnew2=(Xnew2.*(~(FU+FL)))+ub.*FU+lb.*FL;
            CostNew=fobj(Xnew2);            
            
            if CostNew<cost(i)
                X(i,:)=Xnew2;
                cost(i,:)=CostNew;
            else
                if rand<w(randi(dim)) 
                    SM = RungeKutta(X(i,:),Xnew2,DelX);
                    Xnew = (Xnew2-rand.*Xnew2)+ SF(i)*(SM+(2*rand(1,dim).*Best_X-Xnew2));  % (Eq. 20)
                    
                    FU=Xnew>ub;FL=Xnew<lb;Xnew=(Xnew.*(~(FU+FL)))+ub.*FU+lb.*FL;
                    CostNew=fobj(Xnew);
                    
                    if CostNew<cost(i)
                        X(i,:)=Xnew;
                        cost(i,:)=CostNew;
                    end
                end
            end
        end
% End of ESQ 

%% Determine the Best Solution
        if cost(i)<Best_Cost
            Best_X=X(i,:);
            Best_Cost=cost(i);
        end
    end
    
% Save Best Solution at each iteration
Convergence_curve(:,it) = Best_Cost;

disp(['it : ' num2str(it) ', Best Cost 1 = ' num2str(Convergence_curve(1,it)) ', Best Cost 2 = ' num2str(Convergence_curve(2,it))]);


end
disp(['The best optimal values are Minimum Noise: ',num2str(Convergence_curve(1,MaxIt)) ' Minimum Distortion: ' num2str(Convergence_curve(2,MaxIt))]);
end


% A funtion to determine a random number 
%with uniform distribution (unifrnd function in Matlab) 
function z=Unifrnd(a,b,c,dim)
a2 = a/2;
b2 = b/3;
mu = a2+b2;
sig = b2-a2;
z = mu + sig .* (2*rand(c,dim)-1);
end

% A function to determine three random indices of solutions
function [A,B,C]=RndX(nP,i)
Qi=randperm(nP);Qi(Qi==i)=[];
A=Qi(1);
B=Qi(2);
C=Qi(3);
end
