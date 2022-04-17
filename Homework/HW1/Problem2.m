clc; clear; close all;
%% Global Variable

% Variable
% number of vertices
N = 21;

%Density
densMetal = 7000;
densFluid = 1000;
dDensity = densMetal - densFluid;

% Rod Properties
RodLength = 0.1; 
RodRadius = 0.001;
deltaL = RodLength / (N-1);

% Spheres Properties
dL= RodLength/(N-1);
R_rest = dL/10;
R_mid = 0.025;

m_rest = 4/3*pi()*densMetal*R_rest^3;
m_mid = 4/3*pi()*densMetal*R_mid^3;

% Time step
dt = 1e-2;  % second

% Young's modulus
Y = 1e9;

% Gravity
g = 9.81;

% viscosity
visc = 1000;

% total Time
totalTime = 10;

% Number of Quantities
ne = N -1; % Number of Edges
EI = Y * pi() * RodRadius^4 / 4;
EA = Y * pi() * RodRadius^2;

% Geometry
nodes = zeros(N,2);
for c = 1:N
    nodes(c,1) = (c-1)*deltaL; 
end

% Mass Matrix
M = zeros(2*N, 2*N);

for i = 1:2:2*N
    M(i,i) = 4/3*pi()*R_rest^3*densMetal;
    M(i+1,i+1) = 4/3*pi()*R_rest^3*densMetal;
end
M(N,N) = 4/3*pi()*R_mid^3*densMetal;
M(N+1,N+1) = 4/3*pi()*R_mid^3*densMetal;

% Viscous damping matrix
C = zeros(2*N,2*N);
C_rest = 6*pi()*visc*R_rest;
C_mid = 6*pi()*visc*R_mid;

for i = 1:2*N
    C(i,i) = C_rest;
end
C(N,N) = C_mid;
C(N+1,N+1) = C_mid;

% Gravity
W = zeros(2*N,1);
for i = 2:2:2*N
    W(i,1) = -4/3*pi()*R_rest^3*dDensity*g;
end
W(N+1)= -4/3*pi()*R_mid^3*dDensity*g;

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

tol = EI / RodLength^2 * 1e-3;
% Time Matching
for c = 2:N_steps
    fprintf('Time = %f\n', (c-1)*dt);
    q = q0;     % Guess
    
    % Newton Raphson
    err = 10 *tol;
    while err > tol
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
        f = f+ C* (q-q0)/dt;
        J = J +C/dt;
        % Weight
        f = f-W;

        % Update
        q = q-J\f;
        err = sum(abs(f));
    
        
    end
    % update
    u = (q -q0)/dt;     %velocity
    q0 = q;
    
    figure(1);
    plot(q(1:2:end), q(2:2:end), 'ro-');
    axis equal
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