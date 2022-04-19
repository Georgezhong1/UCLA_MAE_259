clc; clear; close all;
% Global Variable
N = 21;
dt = 1e-2;  % second
totalTime = 10;
N_steps = round(totalTime / dt);

% Simulation
[all_yq, all_yv] = P2_simulation(N, dt, totalTime); % simulation
figure(1);
for i = 1:1:length(all_yq)
    plot(all_yq(1:2:end,i), all_yq(2:2:end,i), 'ro-');
    axis equal
    grid on;
    drawnow
    
end

%% Q1: Vertical position and velocity of middle node
all_mid_yq = all_yq(N+1,:); % Middle Node position
all_mid_yv = all_yv(N+1,:); % Middle Node velocity

figure (2);
timeArray =  (1:N_steps)*dt;
plot(timeArray, all_mid_yq, 'k-');
xlabel('Time, t [sec]');
ylabel('Vertical Position of Mid-node, p [meter]');
title('Time VS Vertical Position of Mid-node')
grid on;

figure (3);
timeArray =  (1:N_steps)*dt;
plot(timeArray, all_mid_yv, 'k-');
xlabel('Time, t [sec]');
ylabel('Velocity of Mid-Node, v [meter/sec]');
title('Time VS Velocity of Mid-Node')
grid on;

%% Q2: final deformed shape of the beam
figure (4);
plot(all_yq(1:2:end,length(all_yq)), all_yq(2:2:end,length(all_yq)), 'ro-');
xlabel('X (m)');
ylabel('Y (m)');
axis equal
title('Final Deformation')
grid on;

%% Q3: terminal velocity vs. the number of nodes
num_Node = 20;
dt = 1e-2;
all_yq_c_n = cell(num_Node,1);
all_yv_c_n = cell(num_Node,1);

for N = 3:2:num_Node*2+1
    [all_yq_c_n{N-2}, all_yv_c_n{N-2}] = P2_simulation(N, dt, totalTime);
end

%
Nodes_cell = 3:2:num_Node*2+1;
terminal_velo = zeros(num_Node,1);

%%
for i = 1:2:num_Node*2
    [wid,len] = size(all_yq_c_n{i});
    terminal_velo(i) = all_yv_c_n{i}(wid/2+1,len);
end
terminal_velo = nonzeros(terminal_velo);

figure (5);
axis equal;
plot(Nodes_cell,terminal_velo, 'ro-');
xlabel('Number of Nodes');
ylabel('Terminal Velocity [meter/sec]');
title('Terminal Velocity vs. Number of Nodes')
grid on;

%% Q3: terminal velocity vs. Step_size
num_step = 4;
N = 21;
all_yq_c_s = cell(num_step,1);
all_yv_c_s = cell(num_step,1);

for i = 1:num_step
    dt = 1 * 10^(-i);
    [all_yq_c_s{i}, all_yv_c_s{i}] = P2_simulation(N, dt, totalTime);
end

Nodes_cell = 1:num_step;
terminal_velo = zeros(num_step,1);

%%
for i = 1:num_step
    [wid,len] = size(all_yq_c_s{i});
    terminal_velo(i) = all_yv_c_s{i}(wid/2+1,len);
end

figure (6);
plot(Nodes_cell,terminal_velo, 'ro-');
xlabel('Step size (10^(^-^N^))');
ylabel('Terminal Velocity [meter/sec]');
title('Terminal Velocity vs. Step_size')
grid on;
