clc; clear; close all;

%% Globel Variable

%% Initial Parameters
nv = 50;                % # of vertices
dt = 0.01;              % Time Step
RodLength = 0.2;        % rod length
rho = 1000;             % density
r0 = 0.001;             % rod radius cross section
Y = 10^7;               % Young's modulus
nu = 0.5;               % Poisson Ratio
G = Y / (2*(1+nu));     % Shear modulus
natR = 0.02;            % Natural curvature radius 
g = [0;0;-9.81]';        % Gravity

toler = 1e-3;         % Force Function Tolerance
maximum_iter = 100;     % maximum iteration
totalTime = 5;          % Time

% Stiffness properties
EI = Y *pi() * r0^4/4;  % Bending Stiffness
EA = Y *pi() * r0^2;    % Stretching stiffness
GJ = G *pi() * r0^4/2;  % twisting stiffness

% 
ndof = 4*nv - 1;        % # of Degrees  ?why -1
ne = nv -1;             % # of edges
dm = pi() * r0^2*RodLength * rho / nv;
ScaleSolve = EI / RodLength^2; % Charateristic bending

%% Initial Rod Geometry
nodes = zeros(nv, 3);
dTheta = (RodLength / ne) / natR;
mass = zeros(ndof, 1);
gravity = zeros (ndof, 1);
for c = 1:nv
    % Geometry
    nodes(c,1) = natR * cos((c-1) * dTheta);    % X axis
    nodes(c,2) = natR * sin((c-1) * dTheta);    % Y axis value
    
    % Mass
    if c==1 || c == nv
        mass(4*(c-1) + 1 : 4*(c-1)+3) = dm;
    else
        mass(4*(c-1) + 1 : 4*(c-1)+3) = dm;
    end

    % Gravity
    gravity(4*(c-1) + 1: 4*(c-1) +3) = g;
end


for c = 1:ne                    
    mass(4*c) = pi()*r0^2* RodLength * rho / ne * r0^2/2;       % I = 1/2*m*r^2, Add Inertia
    
    % Reference Undeformed Edge length
    dx = nodes(c+1, :) - nodes(c,:);
    refLen = norm(dx);
end
Fg = mass .* gravity;      % Weight

% Undeformed Voronoi Length
voronoiRefLen = zeros(nv, 1);
for c = 1:nv
    if c == 1
        voronoiRefLen(c) = 0.5 * refLen;
    elseif c == nv
        voronoiRefLen(c) = 0.5 * refLen;
    else
        voronoiRefLen(c) = 0.5 * (refLen ...
            + refLen);
    end
end

%% Loop Parameters
mass = diag(mass);
u = zeros(nv-1,3);
v = zeros(nv-1,3);
tangent_old = zeros(nv-1,3);
tangent_new = zeros(nv-1,3);
a1_old = zeros(nv-1,3);
a1_new = zeros(nv-1,3);
a2_old = zeros(nv-1,3);
a2_new = zeros(nv-1,3);
m1 = zeros(nv-1,3);
m2 = zeros(nv-1,3);
qd = zeros(4*nv-1,1);
delta_mk = zeros(nv-2,1);
free_index = 8:4*nv-1;
q_hist = zeros(4*nv-1,totalTime/dt);


q_new = zeros(nv*4-1,1);
for i = 1:nv
    q_new(i*4-3:i*4-1) = nodes(i,:);
end
q_old = q_new;

t1 = (q_new(5:7) - q_new(1:3))/norm((q_new(5:7) - q_new(1:3))); % Direction Vector
u11 = [-t1(2),t1(1),0];

for i = 1:ne
    % Compute Tangent
    dx = nodes(i+1,:) - nodes(i,:);
    tangent_new(i,:) = dx / norm(dx);
    if i == 1
        u(i,:) = [-tangent_new(i,2), tangent_new(i,1),0];
    else
        u(i,:) = parallel_transport(u(i-1,:), ...
            tangent_new(i-1,:), tangent_new(i,:));
    end
    
    % Compute a1, a2
    a1_new(i,:) = u(i,:);
    a2_new(i,:) = cross(tangent_new(i,:), a1_new(i,:));
    
    % Compute v
    v(i,:) = cross(tangent_new(i,:), u(i,:));

    % Compute m1, m2
    thetak = q_new(4*i);
    m1(i,:) =  cos(thetak)*a1_new(i,:) + sin(thetak)*a2_new(i,:);
    m2(i,:) = -sin(thetak)*a1_new(i,:) + cos(thetak)*a2_new(i,:);
end

kappa = zeros((nv-2),2);
for i = 1: (nv-2)
    kappa(i,:) = computekappa(q_new(4*i-3:4*i-1),q_new(4*i+1:4*i+3), ...
        q_new(4*i+5:4*i+7), m1(i,:),m2(i,:), m1(i+1,:),m2(i+1,:));
end

%% Time Loop
for t = 1:(totalTime/dt)       % time Index
    error = 100;
    while error >= toler
        for jj = 1:ne           % joint index
            dx = nodes(jj+1,:) - nodes(jj,:);
            tangent_new(jj,:) = dx / norm(dx);
            tk1 =(q_new(4*jj+1:4*jj+3) - q_new(4*jj-3:4*jj-1))/norm((q_new(4*jj+1:4*jj+3) - q_new(4*jj-3:4*jj-1)));

            
%             if t == 1 || jj ==1
%                 u(jj,:) = [-tangent_new(jj,2), tangent_new(jj,1),0];
%             else
%                     u(i,:) = parallel_transport(u(i-1,:), ...
%                         tangent_new(i-1,:), tangent_new(i,:));
%             end
            if t == 1
                if jj == 1
                    u(jj,:) = u11;
                else
                    u(jj,:) = parallel_transport(u(jj-1,:)', ...
                        tangent_new(jj-1,:)', tk1);
                end
                a1_new(jj,:) = u(jj,:);
            
            else
                if jj == 1
                    u(jj,:)= u11;
                else
                    u(jj,:) = parallel_transport(u(jj-1,:)', ...
                        tangent_new(jj-1,:)', tk1);
                end
                a1_new(jj,:)= parallel_transport(a1_old(jj,:)',tangent_old(jj,:)',tk1);
            end

            a2_new(jj,:) = cross(tk1, a1_new(jj,:));

            % Compute v
            v = cross(tk1, u(jj,:));
            
            % Compute m1, m2
            tangent_new(jj,:)= tk1;

            thetak = q_new(4*jj);
            m1(jj,:) =  cos(thetak)*a1_new(jj,:) + sin(thetak)*a2_new(jj,:);
            m2(jj,:) = -sin(thetak)*a1_new(jj,:) + cos(thetak)*a2_new(jj,:);
        end

        % Compute Reference Twist
            
        for twi = 1:nv-2
            a1k = parallel_transport(a1_new(twi,:)', ...
                tangent_new(twi,:)', tangent_new(twi+1,:)');  % Time parallel transport

            delta_mk(twi) = signedAngle(a1k, a1_new(twi+1,:)', tangent_new(twi,:)');
        end

        [F,Jacob] = Grad_Hs_EE_rod(q_new,EA,EI,GJ,voronoiRefLen,refLen, ...
                kappa,a1_new,a2_new,tangent_new,delta_mk,m1,m2,g,mass,nv);

        F = F+ mass*(((q_new-q_old)/dt-qd)/dt);
        Jini = mass*(1/dt)*(1/dt);
        Jacob = Jacob + Jini;
        F_free = F(free_index);
        J_free = Jacob(free_index,free_index);
        Jinv=inv(J_free);
        deltaq = Jinv * F_free;
        q_new(free_index) = q_new(free_index) - deltaq ;
        error = sum(abs(F_free));
    end

    q_hist(:,t) = q_new;
    qd = (q_new-q_old)*(1/dt);
    q_old = q_new;
    a1_old = a1_new;
    a2_old = a2_new;
    tangent_old = tangent_new;
end

%% Plot
figure()
for i = 1: size(q_hist)
    
    x = q_hist(1:4:199,i);
    y = q_hist(2:4:199,i);
    z = q_hist(3:4:199,i);
    plot3(x,y,z);
    title('Rod Displacement')
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    axis([-0.05,0.05,-0.05,0.05,-0.1,0.05])
    view([0,0])
    drawnow
    
    
end

figure()
plot(q_hist(4*nv-1,:),'om')
title('Displacement History');
xlabel('Time, t[s]')
ylabel('z-coord of last node, [m]')