% Constants

L_x = 1;
L_y = 1;
Re = 100; % Re = 100, 400, 1000, 3200
n_x = 129;
n_y = 129;
dx = L_x/(n_x-1);
dy = L_y/(n_y-1);
n = n_x + 2;

if Re <= 1000
    dt = (dx ^ 2) * (Re/100);

else
    n_x = 257;
    n_y = 257;
    dx = L_x/(n_x-1);
    dy = L_y/(n_y-1);
    n = n_x + 2;
    dt = (dx ^ 2) * (Re/1000);
end    

mu = 1/Re;
u_top = 1;
v_top = 0;
tolerance = 10 ^ -3;

%%

% Gauss-Seidel Function

function [psi_new,omega_new] = gauss_sdl(psi_prev_timestep,omega_prev_timestep, Re_val, n_x_val, n_y_val,dt_val)

    % Constants
    
    L_x = 1;
    L_y = 1;
    Re = Re_val;
    n_x = n_x_val;
    n_y = n_y_val;
    dx = L_x/(n_x-1);
    dy = L_y/(n_y-1);
    dt = dt_val;
    mu = (1/Re) * 1.05;
    u_top = 1;
    alpha_sor = 1.5;

    % Defining a new array
   

    n = size(psi_prev_timestep,1);
    psi_new = psi_prev_timestep;
    omega_new = omega_prev_timestep;
    tolerance = 10 ^ -3;
    psi_sum = 100;
    omega_sum = 100;

   while psi_sum > tolerance || omega_sum > tolerance

       psi_old = psi_new;
       omega_old = omega_new;
       psi_sum = 0;
       omega_sum = 0;

        for i = 2 : n-1
            for j = 2 : n-1

            psi_new(i,j) = (omega_old(i,j) + (1/(dx ^ 2)) * (psi_new(i,j-1) + psi_new(i,j+1)) + (1/(dy ^ 2)) * (psi_new(i-1,j) + psi_new(i+1,j)))/(2/(dx ^ 2) + 2/(dy ^ 2));
            end 
        end

        psi_new = alpha_sor * psi_new + (1 - alpha_sor) * psi_old;

         % Defining the boundary conditions

        psi_new(:,1) = psi_new(:,2);  % Left wall
        
        psi_new(:,n) = psi_new(:,n-1);  % Right wall
        
        psi_new(1,:) = psi_new(2,:);  % Bottom wall
        
        psi_new(n,:) = u_top * dy + psi_new(n-1,:);  % Top wall


        for i = 2: n-1
            for j = 2: n-1

                u = (psi_new(i+1,j) - psi_new(i-1,j))/(2 * dy);
                v = -(psi_new(i,j+1) - psi_new(i,j-1))/(2 * dx);
                dwbydx = (omega_new(i,j+1) - omega_new(i,j-1))/(2 * dx);
                dwbydy = (omega_new(i+1,j) - omega_new(i-1,j))/(2 * dy);
                omega_new(i,j) = omega_new(i,j) + dt * (mu * ((1/(dx^2)) * (omega_new(i,j+1) - 2 * omega_new(i,j) + omega_new(i,j-1)) + (1/(dy^2)) * (omega_new(i+1,j) - 2 * omega_new(i,j) + omega_new(i-1,j)))...
                    - u * (dwbydx) - v * (dwbydy));

            end
        end
        omega_new = alpha_sor * omega_new + (1 - alpha_sor) * omega_old;

            


        for i = 2 : n-1
            for j = 2: n-1

                psi_sum = psi_sum + abs(psi_old(i,j) - psi_new(i,j));
                omega_sum = omega_sum + abs(omega_old(i,j) - omega_new(i,j));

            end
        end 



        % Defining the boundary conditions

        omega_new(:,1) = (-16/(dx ^ 2)) * psi_new(:,2) - omega_new(:,2);  % Left wall
        
        omega_new(:,n) = (-16/(dx ^ 2)) * psi_new(:,n-1) - omega_new(:,n-1);  % Right wall
        
        omega_new(1,:) = (-16/(dy ^ 2)) * psi_new(2,:) - omega_new(2,:);   % Bottom wall
        
        omega_new(n,:) = (-16/(dy ^ 2)) * psi_new(n-1,:) - (8/dy) * u_top - omega_new(n-1,:);  % Top wall


        
   end

end


%%

% Defining the psi and omega matrix

psi = zeros(n,n);
omega = zeros(n,n);


% Defining the boundary conditions

% Left wall

psi(:,1) = psi(:,2);
omega(:,1) = (-16/(dx ^ 2)) * psi(:,2) - omega(:,2);

% Right wall

psi(:,n) = psi(:,n-1);
omega(:,n) = (-16/(dx ^ 2)) * psi(:,n-1) - omega(:,n-1);

% Bottom wall

psi(1,:) = psi(2,:);
omega(1,:) = (-16/(dy ^ 2)) * psi(2,:) - omega(2,:);

% Top wall

psi(n,:) = u_top * dy + psi(n-1,:);
omega(n,:) = (-16/(dy ^ 2)) * psi(n-1,:) - (8/dy) * u_top - omega(n-1,:);


sum_psi = 100;
sum_omega = 100;
it = 0;


%%

% Main loop

tic

while sum_psi > tolerance || sum_omega > tolerance

    it = it + 1;
    psi_copy = psi;
    omega_copy = omega;
    sum_psi = 0;
    sum_omega = 0;

    [psi,omega] = gauss_sdl(psi,omega,Re,n_x,n_y,dt);

    for i = 1:n
        for j = 1:n
            sum_psi = sum_psi + abs(psi_copy(i,j) - psi(i,j));
            sum_omega = sum_omega + abs(omega_copy(i,j) - omega(i,j));
        end    
    end

    % it
     

end   

total_time = dt * it;


elapsedTime = toc; % End timer and get the elapsed time
fprintf('Total time taken: %.2f seconds\n', elapsedTime);

%% Plotting Streamfunction (psi) and Vorticity (omega) with Contour Lines


[x_grid, y_grid] = meshgrid(linspace(0, L_x, n), linspace(0, L_y, n));


figure;
contour(x_grid, y_grid, psi, 200)
colorbar;
title('Streamfunction');
xlabel('x');
ylabel('y');
grid on;
figure;
contour(x_grid, y_grid, omega, 200); 
colorbar;
title('Vorticity');
xlabel('x');
ylabel('y');
grid on;


%%  Postprocessing for Re = 100


if Re == 100

    u = zeros(n, n);
    v = zeros(n, n);
    u(n,:) = 1;
    v(n,:) = 0;


    for i = 2:n-1
        for j = 2:n-1

            u(i, j) = (psi(i+1, j) - psi(i-1, j)) / (2 * dy);
            v(i, j) = -(psi(i, j+1) - psi(i, j-1)) / (2 * dx);
        end
    end


    velocity_magnitude = sqrt(u.^2 + v.^2);


    [x_grid, y_grid] = meshgrid(linspace(0, L_x, n), linspace(0, L_y, n));


    very_low_speed_levels = linspace(0, 0.0005 * max(velocity_magnitude(:)),35);
    low_speed_levels = linspace(0, 0.01 * max(velocity_magnitude(:)),20); 
    high_speed_levels = linspace(0.01 * max(velocity_magnitude(:)), max(velocity_magnitude(:)), 35);
    
    high_speed_levels = high_speed_levels(1:1:end);

    contour_levels = [very_low_speed_levels,low_speed_levels,high_speed_levels];

    figure;
    quiver(x_grid, y_grid, u, v, 'k');
    hold on;
    contour(x_grid, y_grid, velocity_magnitude, contour_levels); 
    hold off;

    title('Streamlines for Re = 100');
    xlabel('x');
    ylabel('y');
    colorbar;
    grid on;
    
    % ----> data for Re = 100;
    
    
    y_values = linspace(0, L_y, n);
    x_values = linspace(0, L_x, n);
    
    
    j_half = round(0.5 / dx) + 1;
    i_half = round(0.5 / dy) + 1;
    
    
    u_profile = zeros(1, n);
    v_profile = zeros(1, n);
    u_profile(n) = 1;
    v_profile(n) = -(psi(i_half, n) - psi(i_half,n-1)) / (dx);
    
    
    
    for i = 2:n-1
        u_profile(i) = (psi(i+1, j_half) - psi(i-1, j_half)) / (2 * dy);
    end
    
   
    for j = 2:n-1
        v_profile(j) = -(psi(i_half, j+1) - psi(i_half, j-1)) / (2 * dx);
    end
    
    
    y_paper_u = [1.0000, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, ...
               0.6172, 0.5000, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, ...
               0.0625, 0.0547, 0.0000];
    
    u_paper = [1.0000, 0.84123, 0.78871, 0.73722, 0.68717, 0.23151, 0.00332, ...
               -0.13641, -0.20581, -0.21090, -0.15662, -0.10150, -0.06434, ...
               -0.04775, -0.04192, -0.03717, 0.00000];
    
    
    x_paper_v = [1.0000, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, ...
                 0.8047, 0.5000, 0.2344, 0.2266, 0.1563, 0.0938, 0.0781, ...
                 0.0703, 0.0625, 0.0000];
    
    v_paper = [0.00000, -0.05906, -0.07391, -0.08864, -0.10313, -0.16914, -0.22445, ...
               -0.24533, 0.05454, 0.17527, 0.17507, 0.16077, 0.12317, ...
               0.10890, 0.10091, 0.09233, 0.00000];
    
   
    figure;
    subplot(2, 1, 1);
    plot(y_values, u_profile, 'b-', 'LineWidth', 1.5); 
    hold on;
    plot(y_paper_u, u_paper, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r'); 
    title('Comparison of u-velocity Profile at x = 0.5 for Re = 100');
    xlabel('y');
    ylabel('u');
    legend('Computed Profile', 'Reference Data (Ghia and Ghia)');
    grid on;
    
   
    subplot(2, 1, 2);
    plot(x_values, v_profile, 'b-', 'LineWidth', 1.5); 
    hold on;
    plot(x_paper_v, v_paper, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r'); 
    title('Comparison of v-velocity Profile at y = 0.5 for Re = 100');
    xlabel('x');
    ylabel('v');
    legend('Computed Profile', 'Reference Data (Ghia and Ghia)');
    grid on;

end

%%  Postprocessing for Re = 400
 

% ------> PLOTTING THE STREAMLINES


if Re == 400


u = zeros(n, n);
v = zeros(n, n);
u(n,:) = 1;
v(n,:) = 0;


for i = 2:n-1
    for j = 2:n-1

        u(i, j) = (psi(i+1, j) - psi(i-1, j)) / (2 * dy);
        v(i, j) = -(psi(i, j+1) - psi(i, j-1)) / (2 * dx);
    end
end


velocity_magnitude = sqrt(u.^2 + v.^2);


[x_grid, y_grid] = meshgrid(linspace(0, L_x, n), linspace(0, L_y, n));


very_low_speed_levels = linspace(0, 0.0005 * max(velocity_magnitude(:)),20);  
low_speed_levels = linspace(0, 0.01 * max(velocity_magnitude(:)),20);
high_speed_levels = linspace(0.01 * max(velocity_magnitude(:)), max(velocity_magnitude(:)), 35);

high_speed_levels = high_speed_levels(1:1:end); 

contour_levels = [very_low_speed_levels,low_speed_levels,high_speed_levels];

figure;
quiver(x_grid, y_grid, u, v, 'k');
hold on;
contour(x_grid, y_grid, velocity_magnitude, contour_levels);
hold off;

title('Streamlines for Re = 400');
xlabel('x');
ylabel('y');
colorbar;
grid on;



% ----> Data for Re = 400

y_values = linspace(0, L_y, n);
x_values = linspace(0, L_x, n);


j_half = round(0.5 / dx) + 1;
i_half = round(0.5 / dy) + 1;


u_profile = zeros(1, n);
v_profile = zeros(1, n);
u_profile(n) = 1;
v_profile(n) = -(psi(i_half, n) - psi(i_half, n-1)) / (dx);

for i = 2:n-1
    u_profile(i) = (psi(i+1, j_half) - psi(i-1, j_half)) / (2 * dy);
end

for j = 2:n-1
    v_profile(j) = -(psi(i_half, j+1) - psi(i_half, j-1)) / (2 * dx);
end

y_paper_u = [1.0000, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, ...
           0.6172, 0.5000, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, ...
           0.0625, 0.0547, 0.0000];

u_paper = [1.0000, 0.75837, 0.68439, 0.61756, 0.55892, 0.29093, 0.16256, ...
           0.02135, -0.11477, -0.17119, -0.32726, -0.24299, -0.14612, ...
           -0.10338, -0.09266, -0.08186, 0.00000];

x_paper_v = [1.0000, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, ...
             0.8047, 0.5000, 0.2344, 0.2266, 0.1563, 0.0938, 0.0781, ...
             0.0703, 0.0625, 0.0000];

v_paper = [0.00000, -0.12146, -0.15663, -0.19254, -0.22847, -0.23827, ...
           -0.44993, -0.38598, 0.05188, 0.30174, 0.30203, 0.28124, ...
           0.22965, 0.20920, 0.19713, 0.18360, 0.00000];

figure;
subplot(2, 1, 1);
plot(y_values, u_profile, 'b-', 'LineWidth', 1.5);
hold on;
plot(y_paper_u, u_paper, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
title('Comparison of u-velocity Profile at x = 0.5 for Re = 400');
xlabel('y');
ylabel('u');
legend('Computed Profile', 'Reference Data (Ghia and Ghia)');
grid on;

subplot(2, 1, 2);
plot(x_values, v_profile, 'b-', 'LineWidth', 1.5);
hold on;
plot(x_paper_v, v_paper, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
title('Comparison of v-velocity Profile at y = 0.5 for Re = 400');
xlabel('x');
ylabel('v');
legend('Computed Profile', 'Reference Data (Ghia and Ghia)');
grid on;


end


%%  Postprocessing for Re = 1000


if Re == 1000

    u = zeros(n, n);
    v = zeros(n, n);
    u(n,:) = 1;
    v(n,:) = 0;


    for i = 2:n-1
        for j = 2:n-1

            u(i, j) = (psi(i+1, j) - psi(i-1, j)) / (2 * dy);
            v(i, j) = -(psi(i, j+1) - psi(i, j-1)) / (2 * dx);
        end
    end


    velocity_magnitude = sqrt(u.^2 + v.^2);


    [x_grid, y_grid] = meshgrid(linspace(0, L_x, n), linspace(0, L_y, n));

    very_low_speed_levels = linspace(0, 0.005 * max(velocity_magnitude(:)),20);
    low_speed_levels = linspace( 0.01 * max(velocity_magnitude(:)), 0.05 * max(velocity_magnitude(:)),20);
    high_speed_levels = linspace(0.01 * max(velocity_magnitude(:)), max(velocity_magnitude(:)), 35); 
    high_speed_levels = high_speed_levels(1:1:end); 

    contour_levels = [very_low_speed_levels,low_speed_levels,high_speed_levels];

   
    figure;
    quiver(x_grid, y_grid, u, v, 'k'); 
    hold on;
    contour(x_grid, y_grid, velocity_magnitude, contour_levels);
    hold off;

    title('Streamlines for Re = 1000');
    xlabel('x');
    ylabel('y');
    colorbar;
    grid on;
    
    % ----> Data for Re = 1000
    
    
    y_values = linspace(0, L_y, n);
    x_values = linspace(0, L_x, n);
    
   
    j_half = round(0.5 / dx) + 1;
    i_half = round(0.5 / dy) + 1;
    
    
    u_profile = zeros(1, n);
    v_profile = zeros(1, n);
    u_profile(n) = 1;
    v_profile(n) = -(psi(i_half, n) - psi(i_half, n-1)) / (dx);
    
    
    for i = 2:n-1
        u_profile(i) = (psi(i+1, j_half) - psi(i-1, j_half)) / (2 * dy);
    end
    
    
    for j = 2:n-1
        v_profile(j) = -(psi(i_half, j+1) - psi(i_half, j-1)) / (2 * dx);
    end
    
    
    y_paper_u = [1.0000, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, ...
                 0.6172, 0.5000, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, ...
                 0.0625, 0.0547, 0.0000];
    
    u_paper = [1.0000, 0.65928, 0.57492, 0.51117, 0.46604, 0.33304, 0.18719, ...
               0.05702, -0.06080, -0.10648, -0.27805, -0.24299, -0.14612, ...
               -0.10338, -0.09266, -0.08186, 0.00000];
    
    
    x_paper_v = [1.0000, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, ...
                 0.8047, 0.5000, 0.2344, 0.2266, 0.1563, 0.0938, 0.0781, ...
                 0.0703, 0.0625, 0.0000];
    
    v_paper = [0.00000, -0.21388, -0.27669, -0.33714, -0.39188, -0.51550, ...
               -0.42665, -0.31966, 0.02526, 0.32235, 0.33075, 0.37095, ...
               0.32627, 0.30353, 0.29012, 0.27485, 0.00000];
    
    
    figure;
    subplot(2, 1, 1);
    plot(y_values, u_profile, 'b-', 'LineWidth', 1.5);
    hold on;
    plot(y_paper_u, u_paper, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
    title('Comparison of u-velocity Profile at x = 0.5 for Re = 1000');
    xlabel('y');
    ylabel('u');
    legend('Computed Profile', 'Reference Data (Ghia and Ghia)');
    grid on;
    
    subplot(2, 1, 2);
    plot(x_values, v_profile, 'b-', 'LineWidth', 1.5);
    hold on;
    plot(x_paper_v, v_paper, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
    title('Comparison of v-velocity Profile at y = 0.5 for Re = 1000');
    xlabel('x');
    ylabel('v');
    legend('Computed Profile', 'Reference Data (Ghia and Ghia)');
    grid on;

end

%% Post processing for Re = 3200


if Re == 3200

    u = zeros(n, n);
    v = zeros(n, n);
    u(n,:) = 1;
    v(n,:) = 0;


    for i = 2:n-1
        for j = 2:n-1

            u(i, j) = (psi(i+1, j) - psi(i-1, j)) / (2 * dy);
            v(i, j) = -(psi(i, j+1) - psi(i, j-1)) / (2 * dx);
        end
    end


    velocity_magnitude = sqrt(u.^2 + v.^2);


    [x_grid, y_grid] = meshgrid(linspace(0, L_x, n), linspace(0, L_y, n));
    very_low_speed_levels = linspace(0, 0.0005 * max(velocity_magnitude(:)),30);
    low_speed_levels = linspace(0, 0.1 * max(velocity_magnitude(:)),30); 
    high_speed_levels = linspace(0.01 * max(velocity_magnitude(:)), max(velocity_magnitude(:)), 35);
    high_speed_levels = high_speed_levels(1:1:end); 

   
    contour_levels = [very_low_speed_levels,low_speed_levels,high_speed_levels];
    figure;
    quiver(x_grid, y_grid, u, v, 'k');
    hold on;
    contour(x_grid, y_grid, velocity_magnitude, contour_levels);
    hold off;

    title('Streamlines for Re = 3200');
    xlabel('x');
    ylabel('y');
    colorbar;
    grid on;
   

    y_values = linspace(0, L_y, n);
    x_values = linspace(0, L_x, n);
    
    
    j_half = round(0.5 / dx) + 1;
    i_half = round(0.5 / dy) + 1;
    
    
    u_profile = zeros(1, n);
    v_profile = zeros(1, n);
    u_profile(n) = 1;
    v_profile(n) = -(psi(i_half, n) - psi(i_half, n-1)) / (dx);
    
    
    for i = 2:n-1
        u_profile(i) = (psi(i+1, j_half) - psi(i-1, j_half)) / (2 * dy);
    end
    
    
    for j = 2:n-1
        v_profile(j) = -(psi(i_half, j+1) - psi(i_half, j-1)) / (2 * dx);
    end
    
    
    y_paper_u = [1.0000, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, ...
                 0.6172, 0.5000, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, ...
                 0.0625, 0.0547, 0.0000];
    
    u_paper = [1.0000, 0.53236, 0.48296, 0.46547, 0.46101, 0.34682, 0.19791, ...
               0.07156, -0.04272, -0.86636, -0.24427, -0.34323, -0.41933, ...
               -0.37827, -0.35344, -0.32407, 0.00000];
    
    
    x_paper_v = [1.0000, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, ...
                 0.8047, 0.5000, 0.2344, 0.2266, 0.1563, 0.0938, 0.0781, ...
                 0.0703, 0.0625, 0.0000];
    
    v_paper = [0.00000, -0.39017, -0.47425, -0.52357, -0.54053, -0.44307, ...
               -0.37401, -0.31184, 0.00999, 0.28188, 0.29030, 0.37119, ...
               0.42768, 0.41906, 0.40917, 0.39560, 0.00000];
    
    
    figure;
    subplot(2, 1, 1);
    plot(y_values, u_profile, 'b-', 'LineWidth', 1.5);
    hold on;
    plot(y_paper_u, u_paper, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
    title('Comparison of u-velocity Profile at x = 0.5 for Re = 3200');
    xlabel('y');
    ylabel('u');
    legend('Computed Profile', 'Reference Data (Ghia and Ghia)');
    grid on;
    
    
    subplot(2, 1, 2);
    plot(x_values, v_profile, 'b-', 'LineWidth', 1.5); 
    hold on;
    plot(x_paper_v, v_paper, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
    title('Comparison of v-velocity Profile at y = 0.5 for Re = 3200');
    xlabel('x');
    ylabel('v');
    legend('Computed Profile', 'Reference Data (Ghia and Ghia)');
    grid on;

end


%%