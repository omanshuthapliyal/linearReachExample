function [lambdas, xStars] = reachSet(A,B,uVerts,x0,lambda0,t,t0)
% Inputs:
% A, B system matrices
% x0 is [2, nh] sizes matrix with each column representing the i-th
% point of contact for the i-th hyperplane

nu = size(uVerts,1); % Dimension of u-vector
nu_planes = size(uVerts,2);    % number of u-hyperplanes
nx = size(lambda0,1); % Dimension of x-vector
nx_planes = size(lambda0,2);   % number of x-hyperplanes
lambdas = zeros(size(lambda0));

% Propagating planes using:
% https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.52.1682&rep=rep1&type=pdf

% Equation 13
for i = 1:nx_planes
    % Generate hyperplane normals for each hyperplane
    lambdas(:,i) = expm(-(t-t0)*A')*lambda0(:,i);
end

% Equation 14
% Finding Optimal Control
uStar = zeros(nu,nx_planes);
for i = 1:nx_planes
    fval = lambdas(:,i)' * uVerts;
    uStar(:,i) = uVerts(:,fval==max(fval));
end

        
% Equation 12
% Propagating points of plane contact
xStars = zeros(size(x0));
t = linspace(t0,t,1e1);
for i = 1:nx_planes
    sys = ss(A,B, eye(nx),0);
    u = repmat(uStar(:,i), 1,numel(t));
    x = lsim(sys, u, t, x0(:,i));
    xStars(:,i) = x(end,:);
end

end

