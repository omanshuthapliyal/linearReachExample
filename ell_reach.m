function [q, Qs, l] = ell_reach(tspan, A, B, q0, Q0, l0, p, P)

% Implemented for continuous time invariant system and input set.
% Input:
% A, B are system and input matrix of the linear system, respectively.
% q0 and Q0 are the center and scale of the initial set at time 0,
% respectively.
% l0 is the initial support vector at time 0.
% p and P are the center and scale of the input set, respectively,
%
% Output: 
% q and Q are the center and scale of the reach set at time
% t, respectively.
% l is the support vector at time tf.

n = size(B, 1);
y0 = [q0; l0; reshape(Q0, n^2, 1); reshape(Q0, n^2, 1)];
[t, y] = ode45(@(t, y) dyn(t, y, A, B, l0, p, P), tspan, y0);
q = reshape(y(:,1:n).', n, []);
l = reshape(y(:,1+n:2*n).', n, []);
Qs = reshape(y(:,1+(n^2+2*n):(2*n^2+2*n)).', n, n, []);


function y_dot = dyn(t, y, A, B, l0, p, P)

n = size(B, 1);

q = reshape(y(1:n), n, 1);
l = reshape(y(1+n:2*n), n, 1);
Q = reshape(y(1+2*n:(n^2+2*n)), n, n);
Qs = reshape(y(1+(n^2+2*n):(2*n^2+2*n)), n, n);

q_dot = A*q + B*p;
l_dot = -A.'*l;

pi = (l.'*B*P*B.'*l)^0.5*(l0.'*Q*l0)^(-0.5);
Q_dot = pi*Q + 1/pi*expm(-A*t)*B*P*B.'*expm(-A.'*t);

Qs_dot = A*Qs + Qs*A.' + pi*Qs + 1/pi*B*P*B.';


y_dot = [q_dot; l_dot; reshape(Q_dot, n^2, 1); reshape(Qs_dot, n^2, 1)];