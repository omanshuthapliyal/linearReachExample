clearvars;
warning off;
close all;
clc;
% A = [0.2,0.15;0,0.5];
B = [0.1,0.1;0,1];
tSim = 2.4;



tt = 42;
A = [cosd(tt) -sind(tt); sind(tt) cosd(tt)];
x0 = [4.5;4.75];
rx = 0.5;

x01 = [-5;-5];
tt = -45;
A1 = [cosd(tt) -sind(tt); sind(tt) cosd(tt)];
B1 = 2*B;


x02 = [-4.5;-1.5];
tt = -45;
A2 = [cosd(tt) -sind(tt); sind(tt) cosd(tt)];
A2 = 1.2*A1;
B2 = 5*B;

% x03 = [-2.5;3.75];
% tt = 30;
% A3 = [cosd(tt) -sind(tt); sind(tt) cosd(tt)];
% B3 = 1.25*B;


% Initial set is a circle of radius rx around x0
% Initial hyperplanes of polytopes are:
% < c0i, x > \leq \gamma0i*
nx_planes = 12;
thetas = linspace(0,2*pi, nx_planes+1);
Gammas = rx.*ones(1,nx_planes);
% Initializing Costates for each polytope plane
Lambda0 = [cos(thetas);sin(thetas)];
X0 = Lambda0.*rx + x0;

Gammas1 = rx.*ones(1,nx_planes);
% Initializing Costates for each polytope plane
Lambda01 = [cos(thetas);sin(thetas)];
X01 = Lambda01.*rx + x01;

Gammas2 = rx.*ones(1,nx_planes);
% Initializing Costates for each polytope plane
Lambda02 = [cos(thetas);sin(thetas)];
X02 = Lambda02.*rx + x02;

%% admissible control is ||u|| <= rho
rho = 0.1;
n_uverts = 12;
thetas = linspace(0,2*pi, n_uverts+1);
uVertices = zeros(2,n_uverts);
rOut = rho/(cos(pi/n_uverts));
for i = 1:n_uverts
    thet = thetas(i);
    uVertices(:,i) = rOut*[cos(thet); sin(thet)];
    %     plot(uVertices(1,i),uVertices(2,i),'b+'), hold on;
end

% Simulating Reach Sets
DT = 0.03;
tVec = [DT:DT:tSim];
Lambda = zeros(size(Lambda0));
Lambda1 = zeros(size(Lambda01));

% polyA = Lambda0';
% polyb = [gamma01;gamma02;gamma03;gamma04];

% X0 = (polyA.*polyb)' + x0;
% V = con2vert(polyA,polyb);
% [k,~] = convhull(V);
% plot(V(k,1),V(k,2),'b-', 'LineWidth', 2), hold on;
subplot(2,2,[1,3]);

[k,~] = convhull(X0');
plot(X0(1,k),X0(2,k),'g-', 'LineWidth', 2); hold on;

[k,~] = convhull(X01');
plot(X01(1,k),X01(2,k),'g-', 'LineWidth', 2); hold on;

[k,~] = convhull(X02');
plot(X02(1,k),X02(2,k),'g-', 'LineWidth', 2); hold on;

xlabel('$x_1$','Interpreter','Latex','FontSize',20)
ylabel('$x_2$','Interpreter','Latex','FontSize',20)
grid on; grid minor;
axis equal
xlim([-30 10]);
ylim([-10 40])
t_prev = 0;

h = [];
for t = 1:numel(tVec)
    t_curr = tVec(t);
    for i = 1:numel(h)
        set(h(i), 'Visible', 'off');
    end
    
    % Plant A, B
    subplot(2,2,[1,3]);
    
    [lambdas, xStars] = reachSet(A,B,uVertices,X0,Lambda0,t_curr,t_prev);
    Lambda = lambdas;
    xx = xStars';
    [k,~] = convhull(xx);
    plot(xx(k,1),xx(k,2),'m-', 'LineWidth', 0.5);
    Lambda0 = Lambda;
    X0 = xStars;
    
    % Plant A1, B1
    [lambdas1, xStars1] = reachSet(A1,B1,uVertices,X01,Lambda01,t_curr,t_prev);
    Lambda1 = lambdas1;
    xx1 = xStars1';
    [k,~] = convhull(xx1);
    plot(xx1(k,1),xx1(k,2),'b-', 'LineWidth', 0.5);
    Lambda01 = Lambda1;
    X01 = xStars1;
    
    % Plant A2, B2
    [lambdas2, xStars2] = reachSet(A2,B2,uVertices,X02,Lambda02,t_curr,t_prev);
    Lambda2 = lambdas2;
    xx2 = xStars2';
    [k,~] = convhull(xx2);
    plot(xx2(k,1),xx2(k,2),'r-', 'LineWidth', 0.5);
    Lambda01 = Lambda1;
    X02 = xStars2;
    
    t_prev = t_curr;
    
    % Intersecting Polygons
    poly = polyshape(xx(k,1),xx(k,2));
    poly1 = polyshape(xx1(k,1),xx1(k,2));
    poly2 = polyshape(xx2(k,1),xx2(k,2));
    
    polySect = intersect(poly,poly1);
    polySect1 = intersect(poly,poly2);
    polySect2 = intersect(poly2,poly1);
    
    [X1(t),Y1(t)] = poly.centroid;
    [X2(t),Y2(t)] = poly1.centroid;
    [X3(t),Y3(t)] = poly2.centroid;
    
    subplot(2,2,[1,3]);
    
    plot(polySect, 'Facecolor', '#7E2F8E', 'FaceAlpha',0.35)
    plot(polySect1, 'Facecolor', '#7E2F8E', 'FaceAlpha',0.35)
    plot(polySect2, 'Facecolor', '#7E2F8E', 'FaceAlpha',0.35)
    b = [1;1;1]; % 2 => unsafe, 1 => safe
    if ~isempty(polySect.Vertices)
        b(1) = 2;
        b(2) = 2;
    end
    if ~isempty(polySect1.Vertices)
        b(1) = 2;
        b(3) = 2;
    end
    if ~isempty(polySect2.Vertices)
        b(3) = 2;
        b(2) = 2;
    end
    
    subplot(2,2,[1,3]);
    
    h(1) = plot([X1(t), X2(t)], [Y1(t),Y2(t)], 'k-', 'LineWidth', 1.5);
    h(2) = plot([X1(t), X3(t)], [Y1(t),Y3(t)], 'k-', 'LineWidth', 1.5);
    h(3) = plot([X3(t), X2(t)], [Y3(t),Y2(t)], 'k-', 'LineWidth', 1.5);
    h(4) = text(X1(t), Y1(t), '$1$','Interpreter','Latex','FontSize',16);
    h(5) = text(X2(t), Y2(t), '$2$','Interpreter','Latex','FontSize',16);
    h(6) = text(X3(t), Y3(t), '$3$','Interpreter','Latex','FontSize',16);
    
    subplot(2,2,4);
    cats = categorical({'Agent 1','Agent 2','Agent 3'});
    
    

    bar(cats, b);
    ylim([0.8 2.2])
    yticks([1 2])
    yticklabels({'Safe','Unsafe'})

    F(t) = getframe(gcf) ;
    pause(0.01)
    
end
% v = VideoWriter('2-agent.avi');
% v.Quality = 95;
videoWriter(F, 10, '2-agent.avi')


