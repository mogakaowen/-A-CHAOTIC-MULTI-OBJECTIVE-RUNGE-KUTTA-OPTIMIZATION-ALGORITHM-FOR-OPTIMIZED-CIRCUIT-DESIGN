function [lb,ub,dim,fobj] = BenchmarkFunctions2(F)

D=30;
switch F
    case 'F1'
        fobj = @F1;
        lb=-100;
        ub=100;
        dim=D;
        
    case 'F2'
        fobj = @F2;
        lb=-100;
        ub=100;
        dim=D;
        
    case 'F3'
        fobj = @F3;
        lb=-4;
        ub=4;
        dim=D;
        
    case 'F4'
        fobj = @F4;
        lb=-3;
        ub=3;
        dim=D;
        
    case 'F5'
        fobj = @F5;
        lb=0;
        ub=pi;
        dim=D;
        
    case 'F6'
        fobj = @F6;
        lb=-100;
        ub=100;
        dim=D;
        
    case 'F7'
        fobj = @F7;
        lb=-100;
        ub=100;
        dim=D;
        
   case 'F8'
        fobj = @F8;
        lb=0;
        ub=1;
        dim=D;
        
    case 'F9'
        fobj = @F9;
        lb=0;
        ub=1;
        dim=D;
        
    case 'F10'
        fobj = @F10;
        lb=0;
        ub=1;
        dim=D;
    case 'F11'
        fobj = @F11;
        lb=-100;
        ub=100;
        dim=D;
        
    case 'F12'
        fobj = @F12;
        lb=-100;
        ub=100;
        dim=D;
        
    case 'F13'
        fobj = @F13;
        lb=-100;
        ub=100;
        dim=D;  
 
    case 'F14'
        fobj = @F14;
        lb=-5;
        ub=5;
        dim=D;
        
    case 'F15'
        fobj = @F15;
        lb=-4;
        ub=4;
        dim=D;
         
end
end


%% Schaffer function
function f = F1(x)
f(1) = x(1)^2;
f(2) = (x(2) - 2)^2;
end

%% Poloni function
function f = F2(x)

f(1) = (1 + (x(1) + 1)^2*(3*x(2)^2 - 14*x(2) + 19) - (x(2) + 1)^2*(3*x(1)^2 - 14*x(1) + 19));
f(2) = (30 + (x(1) - 2)^2 + (x(2) - 2)^2);

end

%% Viennet2 function
function f = F3(x)

x1 = x(1);
x2 = x(2);

f(1) = (x1-2)^2/2 + (x2+1)^2/13 + 3;
f(2) = (x1+x2-3)^2/36 + (-x1+x2+2)^2/8 - 17;

end

%% Viennet3 function
function f = F4(x)

x1 = x(1);
x2 = x(2);

f(1) = 0.5*(x1^2+x2^2) + sin(x1^2+x2^2);
f(2) = (3*x1-2*x2+4)^2/8 + (x1-x2+1)^2/27 + 15;

end

%% Tanaka Test function
function f = F5(x)

x1 = x(1);
x2 = x(2);
term1 = 1 - cos(x1)*cos(x2);
term2 = abs(term1)^0.98;
f(1) = term1 + 0.1*term2;
f(2) = sin(x1) + 0.1*sin(x2);
end

%% ZDT2 function
function f = F6(x)
n = length(x); % number of decision variables

% Objective 1
f(1) = x(1);

% Objective 2
g = 1 + 9 / (n-1) * sum(x(2:end));
f(2) = g * (1 - (x(1) / g)^2);
end

%% ZDT3 function
function f = F7(x)
n = length(x); % number of decision variables

% Objective 1
f(1) = x(1);

% Objective 2
g = 1 + 9 / (n-1) * sum(x(2:end));
f(2) = g * (1 - sqrt(x(1) / g) - (x(1) / g) * sin(10 * pi * x(1)));
end

%% ZDT4 function 
function f = F8(x)
n = length(x);
f(1) = x(1);
g = 1 + 10*(n-1) + sum(x(2:end).^2 - 10*cos(4*pi*x(2:end)));
h = 1 - sqrt(f(1)/g);
f(2) = g*h;
end

%% ZDT5 function 
function f = F9(x)

n = numel(x);

f(1) = 1 + sum(x(2:end));

g = 1 + 9/(n-1)*sum(x(2:end));

h = 1 - sqrt(f(1)/g);

f(2) = g*h;

end

%% ZDT6 function
function f = F10(x)

n = numel(x);

f(1) = 1 - exp(-4*x(1))*sin(6*pi*x(1))^6;

g = 1 + 9*(sum(x(2:end))/(n-1))^0.25;

h = 1 - (f(1)/g)^2;

f(2) = g*h;

end

%% DTLZ1 function
function f = F11(x)

M = 2; % number of objectives
n = length(x); % number of decision variables
k = n - M + 1; % number of equality constraints

% Compute g
g = sum((x(M:end) - 0.5).^2 - cos(20 * pi * (x(M:end) - 0.5)));

% Compute f1
f(1) = 0.5 * prod(x(1:M-1)) * (1 + g);

% Compute f2
f(2) = 0.5 * (1 - x(1)) * (1 + g);

end

%% DTLZ2 function
function f= F12(x)
M = 2; % number of objectives
n = length(x); % number of decision variables
k = n - M + 1; % number of equality constraints

% Compute g
g = sum((x(M:end) - 0.5).^2);

% Compute f1
f(1) = (1 + g) * cos(x(1) * pi / 2) * cos(x(2) * pi / 2);

% Compute f2
f(2) = (1 + g) * cos(x(1) * pi / 2) * sin(x(2) * pi / 2);
end

%% DTLZ4 function
function f = F13(x)
M = 2; % number of objectives
n = length(x); % number of decision variables
k = n - M + 1; % number of equality constraints

% Compute g
g = sum((x(M:end) - 0.5).^2);

% Compute theta
theta = zeros(1, M-1);
for i = 1:M-1
    theta(i) = pi / (4 * (1 + g)) * (1 + 2 * g * x(i));
end

% Compute f1
f(1) = (1 + g) * prod(cos(theta));

% Compute f2
f(2) = (1 + g) * prod(cos(theta(1:end-1))) * sin(theta(end));
end

%% Kursawe function
function f = F14(x)

n = 3;
f(1) = 0;
f(2) = 0;
for i = 1:n-1
    xi = x(i);
    xi1 = x(i+1);
    exp1 = exp(-0.2*sqrt(xi^2+xi1^2));
    exp2 = exp(-0.2*sqrt(xi1^2+x(i+2)^2));
    f(1) = f(1) - 10*exp1 - 10*exp2;
    f(2) = f(2) + abs(xi)^0.8 + 5*sin(xi^3);
end
end

%% Fonseca-Fleming function
function f = F15(x)

x1 = x(1);
x2 = x(2);
term1 = 1 - exp(-((x1-1/sqrt(2))^2 + (x2-1/sqrt(2))^2));
term2 = 1 - exp(-((x1+1/sqrt(2))^2 + (x2+1/sqrt(2))^2));
f(1) = term1 + term2;
f(2) = term1 - term2;
end

