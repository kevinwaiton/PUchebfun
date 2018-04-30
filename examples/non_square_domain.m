%We provide a way to build approximations on non square domains with
%PUFunLS. Here you must provide the function, domain, and bounding box.
%Options can be specified after.
DOMAIN = DoubleAstroid();
F = PUFunLS(@(x,y)atan(x.^2+y),DOMAIN,[-0.95 0.95;-0.95 0.95],'tol',1e-8);

x = linspace(-1,1,100)';
y = linspace(-1,1,100)';

[X,Y] = ndgrid(x,y);
 
V = F.ChebRoot.evalfGrid({x y});

ind = DOMAIN.Interior([X(:) Y(:)]);

B = DOMAIN.Boundary(800);

VB = F.ChebRoot.evalf(B);

P = [X(ind) Y(ind)];

%Figure out triangulation with star boundary
TRI = delaunayTriangulation([B;P],[(1:length(B)-1)' (2:length(B))'; length(B) 1]);

%Determine interior triangles
TF = isInterior(TRI);

defaultOpts = {'facecolor', 'flat', 'edgealpha', .5, 'edgecolor', 'none'};

subplot(1,2,1);

trisurf(TRI.ConnectivityList(TF,:),TRI.Points(:,1),TRI.Points(:,2),[VB;V(ind)],defaultOpts{:});

axis tight; axis square;

subplot(1,2,2);
plot(B(:,1),B(:,2),'LineWidth',2,'color','black'); hold on;
show(F.ChebRoot); hold off;
axis tight; axis square;

