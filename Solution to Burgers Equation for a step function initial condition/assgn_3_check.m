%% Implicit
 
T = 10;  % final time
dt = 0.01; 
dx = 0.08;
x0 = 1; % x range
x = -x0:dx:x0;
g = size(x,2)-1;
t = 0:dt:T;
Re = 10;
U = [];  % store u at each time step and each x
u_old = zeros(size(x));   % u at each x at a specific time 

N = size(x,2);  % number of nodes

% Apply BCs
u_old(1) = 1;
u_old(end) = 0;
u = u_old;

% Apply ICs
u_old(1:floor(N/2)) = 1;
u_old(floor(N/2)+1:end) = 0;

U = [u_old];

% % Building Thomas algo matrix
% D = zeros(g-1,1);
% A = ( (1/dt) + 2/((Re)*(dx^2)) )*eye(g-1);  % b
% for j = 1:g-1
%     D(j) = (u_old(j + 1)/dt) + (u_old(j + 1+1)^2 - u_old(j + 1-1)^2)/(4*(dx^2));
%     if(j-1~=0)
%         A(j-1,j) = - (u_old(j-1 + 1)/((2)dx)) - (1/(Re(dx^2))); % a
%     end
%     if(j+1<=g-1)
%         A(j+1,j) = (u_old(j+1 + 1)/((2)*dx)) - (1/(Re*dx^2)); % c
%     end
% end
% D(1) = D(1) - (- u_old(1)/(2*dx) - 1/(Re*dx^2))*u_old(1); 
% D(g-1) = D(g-1) - (u_old(j+1 + 1)/(2*dx) - 1/(Re*dx^2))*u_old(g+1);


% ADI loop
for n = 1:T/dt

    % Building Thomas algo matrix
    D = zeros(g-1,1);
    A = ( (1/dt) + (2/(Re*(dx^2))) )*eye(g-1);  % b
    for j = 1:g-1

        D(j) = u_old(j + 1)/dt + ((u_old(j + 1+1)^2) - (u_old(j + 1-1)^2))/(4*(dx));

        if(j-1~=0)
            A(j,j-1) = - (u_old(j-1 + 1)/(2*dx)) - (1/(Re*dx^2)); % a
        end
        if(j+1<=g-1)
            A(j,j+1) = (u_old(j+1 + 1)/(2*dx)) - (1/(Re*dx^2)); % c
        end
    end
    D(1) = D(1) - (- u_old(1)/(2*dx) - 1/(Re*(dx^2)))*u_old(1); 
    D(g-1) = D(g-1) - ( u_old(g+1)/(2*dx) - 1/(Re*(dx^2)) )*u_old(g+1);

    u(2:g) = (A\D)';
    
    % % Forward Elemination
    % A_ = eye(g-1);
    % A_(1,2) = A(1,2)/A(1,1);
    % for l = 1+1:g-2
    %     A_(l,l+1) = A(l,l+1)/(A(l,l)-A(l,l-1)*A_(l-1,l));
    % end
    % 
    % % Forward Elemination
    % D_ = D;
    % D_(1) = D(1)/A(1,1);
    % for l = 1+1:g-1
    %     D_(l) = ( D(l) - (A(l,l-1)*D_(l-1)) )/( A(l,l) - (A(l,l-1)*A_(l-1,l)) );
    % end
    % 
    % % Back Substitution
    % u(g) = D_(g-1);
    % for m = 1:g-2
    %     u(g-m) = D_(g-1-m) - (A_(g-m-1,g-m)*u(g-m+1));
    % end

    U = [U;u];
    u_old = u;
end

figure;
h = plot(U(1, :), 'LineWidth', 2); % Plot the first row to initialize
ylim([min(U(:)) - 0.5, max(U(:)) + 0.5]); % Adjust Y-axis based on data range
xlabel('Index');
ylabel('Value');
title('Live Plot of Array Values');

% Loop through each row (time step) and update the plot
for t = 2:size(U, 1) % Start from second row since first is already plotted
    h.YData = U(t, :); % Update the YData of the plot
    drawnow;              % Update the figure window
    pause(dt);           % Pause to simulate live plotting (adjust as needed)
end