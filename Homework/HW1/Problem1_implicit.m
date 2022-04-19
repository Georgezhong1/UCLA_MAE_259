clc; clear; close all;
%% Global Variable

% Variable
% number of vertices
N = 3;

%Density
densMetal = 7000;
densFluid = 1000;
dDensity = densMetal - densFluid;

% Spheres Properties
R1 = 0.005;
R2 = 0.025;
R3 = 0.005;
m1 = 4/3*pi()*densMetal*R1^3;
m2 = 4/3*pi()*densMetal*R2^3;
m3 = 4/3*pi()*densMetal*R3^3;

% Time step
dt = 1e-2;  % second

% Rod Properties
RodLength = 0.1; 
RodRadius = 0.001;
deltaL = RodLength / (N-1);

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
M(1,1) = 4/3*pi()*R1^3*densMetal;
M(2,2) = 4/3*pi()*R1^3*densMetal;
M(3,3) = 4/3*pi()*R2^3*densMetal;
M(4,4) = 4/3*pi()*R2^3*densMetal;
M(5,5) = 4/3*pi()*R3^3*densMetal;
M(6,6) = 4/3*pi()*R3^3*densMetal;

% Viscous damping matrix
C = zeros(6,6);
C1 = 6*pi()*visc*R1;
C2 = 6*pi()*visc*R2;
C3 = 6*pi()*visc*R3;
C(1,1) = C1;
C(2,2) = C1;
C(3,3) = C2;
C(4,4) = C2;
C(5,5) = C3;
C(6,6) = C3;

% Gravity
W = zeros(2*N,1);
W(2)= -4/3*pi()*R1^3*dDensity*g;
W(4)= -4/3*pi()*R2^3*dDensity*g;
W(6)= -4/3*pi()*R3^3*dDensity*g;

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

all_mid_yq(1) = q(4);
all_mid_yv(1) = u(4);

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
        xk = q(1);
        yk = q(2);
        xkpl = q(3);
        ykpl = q(4);
        % l_k = deltaL;
        dF = gradEs (xk, yk, xkpl, ykpl, deltaL, EA);
        dJ = hessEs (xk, yk, xkpl, ykpl, deltaL, EA);
        f(1:4) = f(1:4) + dF;
        J(1:4, 1:4) = J(1:4, 1:4) + dJ;
    
        % spring # 2 between 2-3
        xk = q(3);
        yk = q(4);
        xkpl = q(5);
        ykpl = q(6);
        % l_k = deltaL;
        dF = gradEs (xk, yk, xkpl, ykpl, deltaL, EA);
        dJ = hessEs (xk, yk, xkpl, ykpl, deltaL, EA);
        f(3:6) = f(3:6) + dF;
        J(3:6, 3:6) = J(3:6, 3:6) + dJ;
    
        % bending spring b/t node 1, 2, 3
        xkml = q(1);
        ykml = q(2);
        xk = q(3);
        yk = q(4);
        xkpl = q(5);
        ykpl = q(6);
        curvature0 = 0;
        dF = gradEb (xkml, ykml, xk, yk, xkpl, ykpl, curvature0, deltaL, EI);
        dJ = hessEb (xkml, ykml, xk, yk, xkpl, ykpl, curvature0, deltaL, EI);
        f(1:6) = f(1:6)+dF;
        J(1:6, 1:6) = J(1:6, 1:6)+dJ;
    
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
    title('Final Deformation')
    grid on;
    axis equal
    drawnow

    % storing value
    all_mid_yq(c) = q(4);       %middle sphere velo
    all_mid_yv(c) = u(4)
end
%%
figure (2);
timeArray =  (1:N_steps)*dt;
plot(timeArray, all_mid_yv, 'k-');
xlabel('Time, t [sec]');
ylabel('Velocity of mid-node, v [meter/sec]')
title('Time VS Velocity');
grid on;