clc; clear;
%% Global Variable

% Variable
% number of vertices
N = 50;

%Density
densAlum = 2700;

% Force
P = 2000;

% Rod Properties
l = 1; 
d = 0.75;
r = 0.011;
R = 0.013;
deltaL = l / (N-1);

% Spheres Properties
dL= l/(N-1);
% R_Sphere
m_Sphere = pi()*(R^2 - r^2)*l*densAlum / (N-1);

% Inertia
I = pi()/4*(R^4 - r^4);

% Time step
dt = 1e-3;  % second

% Modulus of Elasticity
E = 70e9;

% Gravity
g = 9.81;

% total Time
totalTime = 1;

% Number of Quantities
ne = N -1; % Number of Edges
EI = E*I;
EA = E*(pi()*(R^2-r^2));

% Geometry
nodes = zeros(N,2);
for c = 1:N
    nodes(c,1) = (c-1)*deltaL; 
end

% Mass Matrix
M = zeros(2*N, 2*N);
for i = 1:2:2*N
    M(i,i) = m_Sphere;
    M(i+1,i+1) = m_Sphere;
end

% Force
F = zeros(2*N,1);
F(round(N*3/4)*2) = -P;

% Initial DOF vector
q0 = zeros(2*N,1);
for c=1:N
    q0 (2*c - 1) = nodes(c,1);
    q0 (2*c) = nodes(c,2);
end

% Pos and Velo
q = q0;     % DOF vector
u = (q - q0) / dt;      % Velo Vector

% Time steps
N_steps = round(totalTime / dt);
all_mid_yq = zeros(N_steps, 1);      % y - position
all_mid_yv = zeros(N_steps, 1);      % y - velocity

all_mid_yq(1) = q(N+1);
all_mid_yv(1) = u(N+1);

tol = EI / l * 1e-3;

%% Time Matching
for c = 2:N_steps
    fprintf('Time = %f\n', (c-1)*dt);
    q = q0;     % Guess
    
    % Newton Raphson
    err = 10 *tol;
    while err > 1e-5
        % Inertia
        f = M/dt * ( ((q-q0)/dt) -u );
        J = M /dt^2;
        
        % Elastic Forces/ Protential
        % spring # 1 b/t 1-2
        for k = 1:(N-1)
            xk = q(2*k-1);
            yk = q(2*k);
            xkpl = q(2*k+1);
            ykpl = q(2*k+2);
            % l_k = deltaL;
            dF = gradEs (xk, yk, xkpl, ykpl, deltaL, EA);
            dJ = hessEs (xk, yk, xkpl, ykpl, deltaL, EA);
            f((2*k-1):(2*k+2)) = f(2*k-1:2*k+2)+dF;
            J(2*k-1:2*k+2, 2*k-1:2*k+2) = J(2*k-1:2*k+2, 2*k-1:2*k+2) + dJ;
        end
        
        % bending spring b/t node 1, 2, 3
        for k = 2:(N-1)
            xkml = q(2*k-3);
            ykml = q(2*k-2);
            xk = q(2*k-1);
            yk = q(2*k);
            xkpl = q(2*k+1);
            ykpl = q(2*k+2);
            curvature0 = 0;
            dF = gradEb (xkml, ykml, xk, yk, xkpl, ykpl, curvature0, deltaL, EI);
            dJ = hessEb (xkml, ykml, xk, yk, xkpl, ykpl, curvature0, deltaL, EI);
            f((2*k-3):(2*k+2)) = f((2*k-3):(2*k+2))+dF;
            J((2*k-3):(2*k+2), (2*k-3):(2*k+2)) = J((2*k-3):(2*k+2), (2*k-3):(2*k+2))+dJ;
        end
    
        % Calculate q using N-R
        % Viscous force
%         f = f+ C* (q-q0)/dt;
%         J = J +C/dt;

        % External Force
        f = f - F;

        % Update
        q(3:(2*N-1)) = q(3:(2*N-1))-J(3:(2*N-1),3:(2*N-1))\f(3:(2*N-1));
        err = sum(abs(f(3:(2*N-1))));
    end
    % update
    u = (q -q0)/dt;     %velocity
    q0 = q;
    
    figure(1);
    plot(q(1:2:end), q(2:2:end), 'ro-');
    axis([0 1 -0.3 0])
    drawnow

    % storing value
    all_mid_yq(c) = q(N+1);       %middle sphere velo
    all_mid_yv(c) = u(N+1);
end
%%
figure (2);
timeArray =  (1:N_steps)*dt;
plot(timeArray, all_mid_yv, 'k-');
xlabel('Time, t [sec]');
ylabel('Velocity of mid-node, v [meter/sec]')