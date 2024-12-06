Re = 10; % Reynolds Number = 10 / 50
dx = 0.08;
dt = 0.01;
U = 0.1;
CFL_num = (U * dt) / dx;   % it should be set less than 1.
x_max = 2;
n = (x_max/dx) + 1;
mu = 1/Re;
iterations = 150;


% Exact Solution
u_exactfunc = @(x, t) 1 - 1/(1 + (erfc(x/sqrt(4 * mu * t))/erfc((x)/sqrt(4 * mu * t))) * exp(-0.5/mu *(x - 0.5 * t)));

% u_exact = zeros(n,1);
% 
% 
% for i = 1 : (n-1)/2
%     u_exact(i) = 1;
% end 
% for i = (n+1)/2 : n
%     u_exact(i) = 0;
% end
% 
% 
% u_exact(1) = 1;
% u_exact(n) = 0;
% 
% it = 0;
% 
% 
% x = linspace(-x_max/2, x_max/2, n);
% figure; 
% hold on;
% title('Analytical Solution');
% xlabel('Position (x)');
% ylabel('u_{analytical}');
% plotHandle_analytical = plot(x, u_exact, 'LineWidth', 2); % Initial plot
% ylim([0 1.2]);
% 
% 
% while it < iterations
% 
%     it = it + 1;
%     t = it *dt;
% 
%     for j = 2 : n-1
% 
%         x = (j-1) * dx -x_max/2;
%         u_exact(j) = u_exactfunc(x,t);
%     end    
% 
% 
%     set(plotHandle_analytical, 'YData', u_exact); 
%     pause(0.05); 
% 
% end



% Method - 1  ----> Explicit Scheme vs Exact Solution

u_exact = zeros(n,1);


for i = 1 : (n/2)
    u_exact(i) = 1;
end 
for i = (n/2) + 1: n
    u_exact(i) = 0;
end


u_exact(1) = 1;
u_exact(n) = 0;

u_explicit = zeros(n,1);

for i = 1 : (n/2)
    u_explicit(i) = 1;
end 
for i = (n/2) + 1: n
    u_explicit(i) = 0;
end


u_explicit(1) = 1;
u_explicit(n) = 0;

x = linspace(-x_max/2, x_max/2, n);
figure;
hold on;
title('Explicit Solution vs Analytical Solution');
xlabel('Position (x)');
ylabel('u');
plotHandle_explicit = plot(x, u_explicit,'r', 'LineWidth', 2, 'DisplayName', 'Explicit'); 
plotHandle_exact = plot(x, u_exact,'b','LineWidth', 2, 'DisplayName', 'Exact'); 
legend;
ylim([0 1.2]); 

it = 0;

while it < iterations
    it = it + 1;
    u_copy = u_explicit;
    k_matrix = zeros(n,0);
    t = it *dt;

    for j = 2:n-1
        k_matrix(j) = (4 * (dx/dt)) * u_explicit(j) + u_explicit(j-1) ^ 2 - u_explicit(j+1) ^ 2 + (4/(Re * dx)) * (u_explicit(j-1) - 2 * u_explicit(j) + u_explicit(j+1));
    end

    for j = 2 : n-1
        u_explicit(j) = (dt/(4 * dx)) * k_matrix(j);
    end  

    for j = 2 : n-1

        x = (j-1) * dx - (x_max/2);
        u_exact(j) = u_exactfunc(x,t);
    end

    set(plotHandle_explicit, 'YData', u_explicit); 
    set(plotHandle_exact, 'YData', u_exact); 
    pause(0.05); 

end



% Method - 2  ----> Implicit Scheme vs Exact Solution



u_exact = zeros(n,1);


for i = 1 : (n/2)
    u_exact(i) = 1;
end 
for i = (n/2) + 1: n
    u_exact(i) = 0;
end

u_exact(1) = 1;
u_exact(n) = 0;




u_implicit = zeros(n,1);

for i = 1 : (n/2)
    u_implicit(i) = 1;
end 
for i = (n/2) + 1: n
    u_implicit(i) = 0;
end


u_implicit(1) = 1;
u_implicit(n) = 0;

it = 0;


x = linspace(-x_max/2, x_max/2, n); 
figure; 
hold on;
title('Implicit Solution vs Analytical Solution');
xlabel('Position (x)');
ylabel('u');
plotHandle_implicit = plot(x, u_implicit,'r','LineWidth', 2,'DisplayName','Implicit');
plotHandle_exact = plot(x, u_exact,'b','LineWidth', 2, 'DisplayName', 'Exact');
legend;
ylim([0 1.2]); 


while it < iterations

    it = it + 1; 
    u_copy = u_implicit;
    t = it * dt;

    a_mat = zeros(n-2,1);
    b_mat = zeros(n-2,1);
    c_mat = zeros(n-2,1);
    d_mat = zeros(n-2,1);

    for j = 2:n-1
        b_mat(j-1) = (1/dt) + (2/(Re * (dx)^2));  
        c_mat(j-1) = u_implicit(j+1)/(2*(dx)) - (1/(Re * (dx)^2));
        a_mat(j-1) = -(u_implicit(j-1)/(2*dx) + 1/(Re * (dx)^2));

        if j == 2
            d_mat(j-1) = (u_implicit(j)/dt) + (u_implicit(j+1)^2 - u_implicit(j-1)^2)/(4*dx) - a_mat(j-1);
        elseif j == n-1
            d_mat(j-1) = (u_implicit(j)/dt) + (u_implicit(j+1)^2 - u_implicit(j-1)^2)/(4*dx) + 0;
        else
            d_mat(j-1) = (u_implicit(j)/dt) + (u_implicit(j+1)^2 - u_implicit(j-1)^2)/(4*dx);
        end    
    end 

    c_mat(1) = c_mat(1)/b_mat(1);
    d_mat(1) = d_mat(1)/b_mat(1);

    for k = 2 : n-2
        c_mat(k) = c_mat(k)/(b_mat(k) - a_mat(k) * c_mat(k-1));
        d_mat(k) = (d_mat(k) - a_mat(k) * d_mat(k-1))/(b_mat(k) - a_mat(k) * c_mat(k-1));
    end 

    
    u_implicit(n-1) = d_mat(n-2);

    for j = n-2 : -1 : 2
        u_implicit(j) = d_mat(j-1) - u_implicit(j+1) * c_mat(j-1);
    end    

     for j = 2 : n-1

        x = (j-1) * dx - (x_max/2);
        u_exact(j) = u_exactfunc(x,t);
    end



    set(plotHandle_implicit, 'YData', u_implicit);
    set(plotHandle_exact, 'YData', u_exact); 
    pause(0.05);

end




% Method - 3  ----> Crank-Nicolson Scheme vs Exact Solution



u_exact = zeros(n,1);


for i = 1 : (n/2)
    u_exact(i) = 1;
end 
for i = (n/2) + 1: n
    u_exact(i) = 0;
end


u_exact(1) = 1;
u_exact(n) = 0;

 

u_cn = zeros(n,1);


for i = 1 : (n/2)
    u_cn(i) = 1;
end 
for i = (n/2) + 1: n
    u_cn(i) = 0;
end


u_cn(1) = 1;
u_cn(n) = 0;

it = 0;
ct = 0;


x = linspace(-x_max/2, x_max/2, n); 
figure;
hold on;
title('Crank-Nicolson Solution vs Analytical Solution');
xlabel('Position (x)');
ylabel('u');
plotHandle_cn = plot(x, u_cn,'r','LineWidth', 2,'DisplayName','Crank-Nicolson');
plotHandle_exact = plot(x, u_exact,'b','LineWidth', 2, 'DisplayName', 'Exact'); 
legend;
ylim([0 1.2]);

while it < iterations
   
    it = it + 1; 
    u_copy = u_cn;
    t = it * dt;

    a_mat = zeros(n-2,1);
    b_mat = zeros(n-2,1);
    c_mat = zeros(n-2,1);
    d_mat = zeros(n-2,1);

    for j = 2:n-1
        b_mat(j-1) = (1/dt) + (1/(Re * (dx)^2));  
        c_mat(j-1) = u_cn(j+1)/(4*(dx)) - (1/(2 * Re * (dx)^2));
        a_mat(j-1) = -(u_cn(j-1)/(4*dx) + 1/(2 * Re * (dx)^2));

        if j == 2
            d_mat(j-1) = (u_cn(j)/dt) + (1/(2 * Re * (dx^2)))*(u_cn(j+1) - 2 * u_cn(j) + u_cn(j-1)) - a_mat(j-1);
        elseif j == n-1
            d_mat(j-1) = (u_cn(j)/dt) + (1/(2 * Re * (dx^2)))*(u_cn(j+1) - 2 * u_cn(j) + u_cn(j-1));
        else
            d_mat(j-1) = (u_cn(j)/dt) + (1/(2 * Re * (dx^2)))*(u_cn(j+1) - 2 * u_cn(j) + u_cn(j-1));
        end    
    end 

    c_mat(1) = c_mat(1)/b_mat(1);
    d_mat(1) = d_mat(1)/b_mat(1);

    for k = 2 : n-2
        c_mat(k) = c_mat(k)/(b_mat(k) - a_mat(k) * c_mat(k-1));
        d_mat(k) = (d_mat(k) - a_mat(k) * d_mat(k-1))/(b_mat(k) - a_mat(k) * c_mat(k-1));
    end 

    u_cn(n-1) = d_mat(n-2);

    for j = n-2 : -1 : 2
        u_cn(j) = d_mat(j-1) - u_cn(j+1) * c_mat(j-1);
    end    

     for j = 2 : n-1

        x = (j-1) * dx - (x_max/2);
        u_exact(j) = u_exactfunc(x,t);
    end

    set(plotHandle_cn, 'YData', u_cn); 
    set(plotHandle_exact, 'YData', u_exact); 
    pause(0.05);

end


