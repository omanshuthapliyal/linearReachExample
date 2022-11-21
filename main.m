tspan = 0:0.03:2.4;
N_vert = 12; % number of approximate ellipses, each corresponding to a initial support vector

eps = 0.5; % radius of initial set
mu = 0.1; % radius of input set 

tt_sys = [42, -45, -45];
b = [1, 2, 5];
q0 = [[4.5; 4.75], [-5; -5], [-4.5; -1.5]]; % center of initial set

sys(3).B = [];
for i = 1:3
    tt = tt_sys(i);
    sys(i).A = [cosd(tt) -sind(tt); sind(tt) cosd(tt)];
    sys(i).B = b(i)*[0.1, 0.1; 0, 1];
    sys(i).q0 = q0(:, i);
end

sys(3).A = sys(3).A*1.2;
Q0 = eps^2*eye(2);
p = [0; 0];
P = mu^2*eye(2);
l0 = [0.5; 0.5];

colors = ['m', 'b', 'r'];

figure();


for sys_i = 1:3
    q = zeros(N_vert,2,81);
    Qs = zeros(N_vert,2,2,81);
    for l_i = 0:N_vert-1
        theta = 2*pi*(l_i+1)/N_vert;
%         theta = 2*pi/6;
        l0 = [cos(theta);sin(theta)];
        
        [q_temp, Qs_temp, l_temp] = ell_reach(tspan, sys(sys_i).A, sys(sys_i).B, sys(sys_i).q0, Q0, l0, p, P);
        q(l_i+1, :, :) = q_temp;
        Qs(l_i+1, :,:,:) = Qs_temp;
    end
    
    for i = 1:size(q, 3)
        contour = ell(q(:,:,i), Qs(:,:,:,i));
        plot(contour(1,:), contour(2,:), 'color', colors(sys_i)); hold on;
    %     plot3(i*ones(1,100), contour(1,:), contour(2,:), 'color', 'k'); hold on;
    end
end


xlim([-30 10]);
ylim([-10 40]);

function contour = ell(center, scale)

N_l = size(center, 1);
theta = linspace(0, 2*pi, 100);
circ = [cos(theta); sin(theta)];
contour = [cos(theta); sin(theta)];
for i = 1:100
    x = circ(:, i);
    
    min_dist = inf;
    dir = [];
    for l_i = 1:N_l
        sc = squeeze(scale(l_i, :,:));
        new_dir = (x.'*inv(sc)*x)^(-0.5)*x;
        
        if min_dist > norm(new_dir)
            dir = new_dir;
            min_dist = norm(new_dir);
        end
    end
    contour(:, i) = center(1,:).' - dir;
end
end

