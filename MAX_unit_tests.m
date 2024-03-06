%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grant Johnson
% 2/27/2024
% Test the Maxwellian Projection and correction routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

% Build required quantities:


%%% WENO TESTING %%%
grid_struct.moments_type = "Simple_No_Weno_reconst_fv";
grid_struct.Nx = 1;
grid_struct.Nv = 10; %150;
Nx = grid_struct.Nx;
Nv = grid_struct.Nv;
f_IC = zeros(Nx,Nv);
grid_struct.x_max = 1;
grid_struct.x_min = 0;
grid_struct.Lx = grid_struct.x_max - grid_struct.x_min;
grid_struct.x = linspace(grid_struct.x_min,grid_struct.x_max,grid_struct.Nx);

grid_struct.v_max = 4;
grid_struct.v_min = -grid_struct.v_max;
grid_struct.Lv = grid_struct.v_max - grid_struct.v_min;
grid_struct.v = linspace(grid_struct.v_min,grid_struct.v_max,grid_struct.Nv);
grid_struct.dv = grid_struct.v(2) - grid_struct.v(1);

% Save the grid_struct.to an app
app.grid_struct = grid_struct;

% Max testing
app.m0 = 1;
app.kb = 1;

% Intial moments:
n0 = 1;
u0 = 1; %1; 
T0 = 0.3; %1;

% Make a maxwellian distribtuion
for i = 1:grid_struct.Nv
    f_IC(i) = maxwellian(n0,u0,T0,grid_struct.v(i),app);
end

% Fake t for ref
t = [0,1];
ts = [0.25,0.75];

% Fix the maxwellian
fprintf("First correction:\n")
tic
f_fixed_Picard = fix_max(f_IC,n0,u0,T0,app);
toc

fprintf("\n\nSecond correction:\n")
tic
f_fixed_Newton = fix_max_2(f_IC,n0,u0,T0,app);
toc

% Initial projection/fixed
[n1,u1,T1] = moments(f_IC,app);
[n2,u2,T2] = moments(f_fixed_Newton,app);
n_corr = [n1,n2];
u_corr = [u1,u2];
T_corr = [T1,T2];


% Plot the exact solution
subplot(2,3,1)
plot(grid_struct.v,f_IC,"-*blue",LineWidth=2)
hold on
plot(grid_struct.v,f_fixed_Newton,":",LineWidth=2)
title("f(v)")
xlabel("v")
ylabel("f(v)")
legend("IC","Fixed")


% Plot the moments
subplot(2,3,2)
plot(t,n0 + t*0,"black")
hold on
plot(ts,n_corr,"*")
title("n")
xlabel("Left: IC, Right Fixed")
ylabel("n")
legend("Specified","Fixed")

% Plot the moments
subplot(2,3,3)
plot(t,u0 + t*0,"black")
hold on
plot(ts,u_corr,"*")
title("u")
xlabel("Left: IC, Right Fixed")
ylabel("u")
legend("Specified","Fixed")

% Plot the moments
subplot(2,3,4)
plot(t,T0 + t*0,"black")
hold on
plot(ts,T_corr,"*")
title("T")
xlabel("Left: IC, Right Fixed")
ylabel("T")
legend("Specified","Fixed")

% Plot the moments
subplot(2,3,5)
plot(grid_struct.v,f_IC - f_fixed_Newton,LineWidth=2)
title("f(v) [IC - fixed]")
xlabel("v")
ylabel("f(v)")
grid_struct.on