F1 = PUchebfun(@(x,y) log(1+10^5*((x-0.5).^2+(y-0.5).^2)),[0 1;0 1]);

%Functions can be computed on a list of points [X Y];
%must be columns
X = [0.25;0.25];
Y = [0.75;0.5];

EF = F2.evalf([X Y]);

%Functions can be computed on a tensor grid as well
%must be columns
x = linspace(0,1,200).';
y = linspace(0,1,200).';

EF_G = F2.evalfGrid({x y});

%We can compute the derivatives, divergence and laplacian
F2 = PUchebfun(@(x,y) log(1+500*((x-0.5).^2+(y-0.5).^2)),[0 1;0 1]);
% Laplacian
L1 = lap(F1);

%First derivative in x
F1dx = diff(F1,1,1);

%We can integrate F2 over the domain
I = sum(F2);

%This allows us to compute the L^2 norm of functions
N = norm(F2);

%We can also numerical combine approximations
F_plus = F1+F2;
F_times = F1*F2;