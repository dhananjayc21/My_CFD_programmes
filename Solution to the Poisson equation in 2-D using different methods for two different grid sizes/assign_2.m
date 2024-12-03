
% Functions for different numerical schemes

function [resultvec,result] = gauss_sdl(phi,s_phi,h)
   
    n = size(phi, 1); % No. of points on the x-axis
    phi_old = phi;
    phi_new = phi;
    tolerance = 10 ^ -6;
    iterations_g_sdl = 0;
    residual = 100; % bigger value given to enter the loop
    residuevec = [];

     while tolerance < residual 
       residual = 0;  
       iterations_g_sdl = iterations_g_sdl + 1;

       for i = n-1 : -1 : 2
           for j = 2 : n-1
                
               phi_new(i,j) = (phi_new(i,j-1) + phi_old(i,j+1) + phi_old(i-1,j) + phi_new(i+1,j) - (h*h) * s_phi((j-1) * h, (n-i) * h))/4;

          end 
       end

       % % Updating the boundary conditions ----> in this problem our bc are dirichlet
       % 
       %  for i = 1 : n
       %      phi_new(i,1) = 2*(b_0y((n-i) * h)) - phi_new(i,2); % left bdry ---> x = 0
       %  end    
       %  for i = 1 : n
       %      phi_new(i,n) = 2*(b_1y((n-i) * h)) - phi_new(i,n-1); % right bdry ---> x = 1
       %  end 
       %  for j = 1 : n
       %      phi_new(1,j) = 2*(b_x0((j-1) * h)) - phi_new(2,j); % top bdry ---> y = 1
       %  end 
       %  for j = 1 : n
       %      phi_new(n,j) = 2*(b_x1((j-1) * h)) - phi_new(n-1,j); % bottom bdry ----> y = 0
       %  end 

        % Calculating the residual

       for i = 1 : n
           for j = 1 : n
                residual = residual +  abs(phi_new(i,j) - phi_old(i,j));
           end
       end

       residuevec(end + 1) = residual;
 
    
       phi_old = phi_new;
     end

     display(iterations_g_sdl);
     result = phi_new;
     resultvec = residuevec;
end
function [resultvec,result] = jacobi(phi,s_phi,h)

    n = size(phi, 1); % No. of points on the x-axis
    phi_old = phi;
    phi_new = phi;
    tolerance = 10 ^ -6;
    residual = 100; % bigger value given to enter the loop
    residuevec = [];
    iterations_jacobi = 0;

     while tolerance < residual 
       residual = 0;  
       iterations_jacobi = iterations_jacobi + 1;

       for i = 2 : n-1
           for j = 2 : n-1
                
               phi_new(i,j) = (phi_old(i,j-1) + phi_old(i,j+1) + phi_old(i-1,j) + phi_old(i+1,j) - (h^2) * s_phi((j-1) * h, (n-i) * h))/4;

          end 
       end

       % % Updating the boundary conditions ----> in this problem our bc are dirichlet
       % 
       %  for i = 1 : n
       %      phi_new(i,1) = 2*(b_0y((n-i) * h)) - phi_new(i,2); % left bdry ---> x = 0
       %  end    
       %  for i = 1 : n
       %      phi_new(i,n) = 2*(b_1y((n-i) * h)) - phi_new(i,n-1); % right bdry ---> x = 1
       %  end 
       %  for j = 1 : n
       %      phi_new(1,j) = 2*(b_x1((j-1) * h)) - phi_new(2,j); % top bdry ---> y = 1
       %  end 
       %  for j = 1 : n
       %      phi_new(n,j) = 2*(b_x0((j-1) * h)) - phi_new(n-1,j); % bottom bdry ----> y = 0
       %  end 


        % Calculating the residual

        for i = 2 : n-1
           for j = 2 : n-1
                residual = residual +  abs(phi_new(i,j) - phi_old(i,j));
           end
        end

       residuevec(end + 1) = residual;
    
       phi_old = phi_new;

     end

     display(iterations_jacobi);
     result = phi_new;
     resultvec = residuevec;
     
end
function [resultvec,result] = sor(phi,s_phi,h)

    n = size(phi, 1); % No. of points on the x-axis
    m = size(phi, 2); % No. of points on the y-axis
    phi_old = phi;
    phi_new = phi;
    residual = 100; % bigger value given to enter the loop
    residualvec = [];
    tolerance = 10 ^ -6;
    iterations_sor = 0;
    alpha_sor = 1.75; % set its value between 0 and 2 -----> I get the least iterations at around 1.9.

     while tolerance < residual 
       residual = 0;  
       iterations_sor = iterations_sor + 1;

       for i = n-1 : -1 : 2
           for j = 2 : m-1
                
               phi_new(i,j) = alpha_sor * (phi_new(i,j-1) + phi_old(i,j+1) + phi_old(i-1,j) + phi_new(i+1,j) - (h^2) * s_phi((j-1) * h, (n-i) * h))/4 + (1 - alpha_sor) * phi_old(i,j);

          end 
       end

       % % Updating the boundary conditions ----> in this problem our bc are dirichlet
       % 
       %  for i = 1 : n
       %      phi_new(i,1) = 2*(b_0y((n-i) * h)) - phi_new(i,2); % left bdry ---> x = 0
       %  end    
       %  for i = 1 : n
       %      phi_new(i,n) = 2*(b_1y((n-i) * h)) - phi_new(i,n-1); % right bdry ---> x = 1
       %  end 
       %  for j = 1 : n
       %      phi_new(1,j) = 2*(b_x0((j-1) * h)) - phi_new(2,j); % top bdry ---> y = 1
       %  end 
       %  for j = 1 : n
       %      phi_new(n,j) = 2*(b_x1((j-1) * h)) - phi_new(n-1,j); % bottom bdry ----> y = 0
       %  end 

       % Calculating the residual

       for i = 1 : n
           for j = 1 : m
               residual = residual +  abs(phi_new(i,j) - phi_old(i,j));
           end
       end

       residualvec(end + 1) = residual;

     
       phi_old = phi_new;
     end

     display(iterations_sor);
     result = phi_new;
     resultvec = residualvec;
     
end
function [resultvec,result] = rbpgs(phi,s_phi,h)

    n = size(phi, 1); % No. of points on the x-axis
    m = size(phi, 2); % No. of points on the y-axis
    phi_old = phi;
    phi_new = phi;
    residual = 100; % bigger value given to enter the loop
    residualvec = [];
    tolerance = 10 ^ -6;
    iterations_rbpgs = 0;

     while tolerance < residual 
       residual = 0;  
       iterations_rbpgs = iterations_rbpgs + 1;

       % Updating the red cells using the black cells (old values of black cells)

       for i = n-1 : -2 : 2
           for j = 2 : 2 : m-1
                
               phi_new(i,j) = (phi_old(i,j-1) + phi_old(i,j+1) + phi_old(i-1,j) + phi_old(i+1,j) - (h^2) * s_phi((j-1) * h, (n-i) * h))/4;

          end 
       end

       for i = n-2 : -2 : 2
           for j = 3 : 2 : m-1
                
               phi_new(i,j) = (phi_old(i,j-1) + phi_old(i,j+1) + phi_old(i-1,j) + phi_old(i+1,j) - (h^2) * s_phi((j-1) * h, (n-i) * h))/4;

          end 
       end

       % Updating the black cells using the red cells (new values of red cells)

       for i = n-2 : -2 : 2
           for j = 2 : 2 : m-1
                
               phi_new(i,j) = (phi_new(i,j-1) + phi_new(i,j+1) + phi_new(i-1,j) + phi_new(i+1,j) - (h^2) * s_phi((j-1) * h, (n-i) * h))/4;

          end 
       end

       for i = n-1 : -2 : 2
           for j = 3 : 2 : m-1
                
               phi_new(i,j) = (phi_new(i,j-1) + phi_new(i,j+1) + phi_new(i-1,j) + phi_new(i+1,j) - (h^2) * s_phi((j-1) * h, (n-i) * h))/4;

          end 
       end


       % 
       % % Updating the boundary conditions ----> in this problem our bc are dirichlet
       % 
       %  for i = 1 : n
       %      phi_new(i,1) = 2*(b_0y((i-1) * h)) - phi_new(i,2); % left bdry ---> x = 0
       %  end    
       %  for i = 1 : n
       %      phi_new(i,n) = 2*(b_1y((i-1) * h)) - phi_new(i,n-1); % right bdry ---> x = 1
       %  end 
       %  for j = 1 : n
       %      phi_new(1,j) = 2*(b_x0((j-1) * h)) - phi_new(2,j); % top bdry ---> y = 1
       %  end 
       %  for j = 1 : n
       %      phi_new(n,j) = 2*(b_x1((j-1) * h)) - phi_new(n-1,j); % bottom bdry ----> y = 0
       %  end 

       % Calculating the residual

       for i = 1 : n
           for j = 1 : m
                residual = residual +  abs(phi_new(i,j) - phi_old(i,j));
           end
       end


       residualvec(end + 1) = residual;

       
    
       phi_old = phi_new;
     end

     display(iterations_rbpgs);
     result = phi_new;
     resultvec = residualvec;
     
end
function [resultvec,result] = adi(phi,s_phi,h)

   
    n = size(phi, 1); % No. of points on the x-axis
    phi_old = phi;
    phi_half = phi;
    residual = 100; % bigger value given to enter the loop
    residualvec = [];
    tolerance = 10 ^ -6;
    iterations_adi = 0;

     while tolerance < residual 

       residual = 0;  
       iterations_adi = iterations_adi + 1;


       % ADI STEP - 1 ----> ALONG THE ROW

       for i = n-1 : -1 : 2


           Ai = 1 * ones(n-2,1); % subdiagonal values
           Bi = -4 * ones(n-2,1); % diagonal values
           Ci = 1 * ones(n-2,1); % superdiagonal values
           Di = zeros(n-2,1);


           for j = 2 : n-1

               if j == 2
                    Di(j-1) = (h*h) * s_phi((j-1)*h,(n-i)*h) - phi_half(i+1,j) - phi_half(i-1,j) - phi_half(i,j-1);
               elseif j == (n-1)
                    Di(j-1) = (h*h) * s_phi((j-1)*h,(n-i)*h) - phi_half(i+1,j) - phi_half(i-1,j) - phi_half(i,j+1);
               else       
                    Di(j-1) = (h*h) * s_phi((j-1)*h,(n-i)*h) - phi_half(i+1,j) - phi_half(i-1,j);
               end

           end


           Ci(1) = Ci(1)/Bi(1);
           Di(1) = Di(1)/Bi(1);


           for k = 2 : n-2

               Ci(k) = Ci(k)/(Bi(k) - Ai(k) * Ci(k-1));
               Di(k) = (Di(k) - Ai(k) * Di(k-1))/(Bi(k) - Ai(k) * Ci(k-1));

           end  


           phi_half(i,n-1) = Di(n-2);


           for j = n-2 : -1 : 2

               phi_half(i,j) = Di(j-1) - phi_half(i,j+1) * Ci(j-1);

           end    


       end


       % % Updating the boundary conditions ----> in this problem our bc's are dirichlet so no changes are required
       % 
       %  for i = 1 : n
       %      phi_new(i,1) = 2*(b_0y((i-1) * h)) - phi_old(i,1); % left bdry ---> x = 0
       %  end    
       %  for i = 1 : n
       %      phi_new(i,n) = 2*(b_1y((i-1) * h)) - phi_old(i,n); % right bdry ---> x = 1
       %  end 
       %  for j = 1 : n
       %      phi_new(1,j) = 2*(b_x0((j-1) * h)) - phi_old(1,j); % top bdry ---> y = 1
       %  end 
       %  for j = 1 : n
       %      phi_new(n,j) = 2*(b_x1((j-1) * h)) - phi_old(n,j); % bottom bdry ----> y = 0
       %  end 


        phi_new = phi_half;



        % ADI STEP - 2 ----> ALONG THE COLUMN


        for j = n-1 : -1 : 2

           Ai = 1 * ones(n-2,1); % subdiagonal values
           Bi = -4 * ones(n-2,1); % diagonal values
           Ci = 1 * ones(n-2,1); % superdiagonal values
           Di = zeros(n-2,1);


           for i = 2 : n-1

               if i == 2
                   Di(i-1) = (h*h) * s_phi((j-1)*h,(n-i)*h) - phi_new(i,j+1) - phi_new(i,j-1) - phi_new(i-1,j);
               elseif i == (n-1)
                   Di(i-1) = (h*h) * s_phi((j-1)*h,(n-i)*h) - phi_new(i,j+1) - phi_new(i,j-1) - phi_new(i+1,j);
               else    
                   Di(i-1) = (h*h) * s_phi((j-1)*h,(n-i)*h) - phi_new(i,j+1) - phi_new(i,j-1);
               end

           end

           Ci(1) = Ci(1)/Bi(1);
           Di(1) = Di(1)/Bi(1);

           for k = 2 : n-2

               Ci(k) = Ci(k)/(Bi(k) - Ai(k) * Ci(k-1));
               Di(k) = (Di(k) - Ai(k) * Di(k-1))/(Bi(k) - Ai(k) * Ci(k-1));

           end  

           phi_new(n-1,j) = Di(n-2);

           for i = n-2 : -1 : 2

               phi_new(i,j) = Di(i-1) - phi_new(i+1,j) * Ci(i-1);

           end    

       end

       % 
       % % % Updating the boundary conditions ----> in this problem our bc are dirichlet
       % % 
       % %  for i = 1 : n
       % %      phi_new(i,1) = 2*(b_0y((i-1) * h)) - phi_new(i,2); % left bdry ---> x = 0
       % %  end    
       % %  for i = 1 : n
       % %      phi_new(i,n) = 2*(b_1y((i-1) * h)) - phi_new(i,n-1); % right bdry ---> x = 1
       % %  end 
       % %  for j = 1 : n
       % %      phi_new(1,j) = 2*(b_x0((j-1) * h)) - phi_new(2,j); % top bdry ---> y = 1
       % %  end 
       % %  for j = 1 : n
       % %      phi_new(n,j) = 2*(b_x1((j-1) * h)) - phi_new(n-1,j); % bottom bdry ----> y = 0
       % %  end



       % Calculating the residual

       for i = 1 : n
           for j = 1 : n
                residual = residual +  abs(phi_new(i,j) - phi_old(i,j));
           end
       end

       % residual
       residualvec(end + 1) = residual;
       
    
       phi_old = phi_new;
     end

     display(iterations_adi);
     result = phi_new;
     resultvec = residualvec;
     
end
function [resultvec,result] = exact(phi,phi_analytical,h)
  
    n = size(phi, 1); % No. of points on the x-axis
    phi_new = phi;
    residualvec = [];


  for i = 2:n-1
      for j = 2:n-1

          phi_new(i,j) = phi_analytical((j-1)*h, (n-i)*h);
      end
  end    

  result = phi_new;
  resultvec = residualvec;

end






% Dimensions of the domain
xlen = 1;
ylen = 1;
n = 41; % can be 41 or 81 for our problem
h = xlen/(n-1); % dx = dy = h

% creating the matrix

phi = zeros(n,n);

% equation ---->  d2(phi)/dx2 + d2(phi)/dy2 = s_phi

s_phi = @(x, y) 2 * sinh(10 * (x - 0.5)) + 40 * (x - 0.5) * cosh(10 * (x - 0.5)) + 100 * ((x - 0.5)^2) * sinh(10 * (x - 0.5)) + 2 * sinh(10 * (y - 0.5)) + 40 * (y - 0.5) * cosh(10 * (y - 0.5)) + 100 * ((y - 0.5)^2) * sinh(10 * (y - 0.5)) + 4 * (x^2 + y^2) * exp(2 * x * y);

phi_analytical = @(x, y) ((x - 0.5)^2) * sinh(10 * (x - 0.5)) + ((y - 0.5)^2) * sinh(10 * (y - 0.5)) + exp (2 * x * y); % ---> analytical solution of our differential equation

% Boundary conditions expressed as functions in which values will be put to fill the boundary of the phi-matrix.

b_0y = @(y) 0.25 * sinh(-5) + ((y - 0.5) ^ 2) * sinh(10 * (y - 0.5)) + 1; % ---> Boundary condition at left wall
b_1y = @(y) 0.25 * sinh(5) + ((y - 0.5) ^ 2) * sinh(10 * (y - 0.5)) + exp(2 * y); % ---> Boundary condition at right wall
b_x0 = @(x) 0.25 * sinh(-5) + ((x - 0.5) ^ 2) * sinh(10 * (x - 0.5)) + 1; % ---> Boundary condition at bottom wall
b_x1 = @(x) 0.25 * sinh(5) + ((x - 0.5) ^ 2) * sinh(10 * (x - 0.5)) + exp(2 * x); % ---> Boundary condition at top wall

% filling the boundary values at the initiation

for i = 1 : n
    phi(i,1) = b_0y((n-i) * h); % left bdry ---> x = 0
end    
for i = 1 : n
    phi(i,n) = b_1y((n-i) * h); % right bdry ---> x = 1
end 
for j = 1 : n
    phi(1,j) = b_x1((j-1) * h); % top bdry ---> y = 1
end 
for j = 1 : n
    phi(n,j) = b_x0((j-1) * h); % bottom bdry ----> y = 0
end 


% These operations I am doing to plot these matrices.

[jacobivec,X_jacobi] = jacobi(phi,s_phi,h);
[sdlvec,X_sdl] = gauss_sdl(phi,s_phi,h);
[sorvec,X_sor] = sor(phi,s_phi,h);
[rbpgsvec,X_rbpgs] = rbpgs(phi,s_phi,h);
[adivec,X_adi] = adi(phi,s_phi,h);
[exactvec,X_exact] = exact(phi,phi_analytical,h);

for i = 1 : n/2
    X_jacobi([i (n-i+1)], :) = X_jacobi([(n-i+1) i], :);
    X_sdl([i (n-i+1)], :) = X_sdl([(n-i+1) i], :);
    X_sor([i (n-i+1)], :) = X_sor([(n-i+1) i], :);
    X_rbpgs([i (n-i+1)], :) = X_rbpgs([(n-i+1) i], :);
    X_adi([i (n-i+1)], :) = X_adi([(n-i+1) i], :);
    X_exact([i (n-i+1)], :) = X_exact([(n-i+1) i], :);
end    

maxlen = max([length(jacobivec), length(sdlvec), length(sorvec), length(adivec)]);
jacobimaxvec = zeros(1,maxlen);jacobimaxvec(1:length(jacobivec)) = jacobivec;
sdlmaxvec = zeros(1,maxlen);sdlmaxvec(1:length(sdlvec)) = sdlvec;
sormaxvec = zeros(1,maxlen);sormaxvec(1:length(sorvec)) = sorvec;
rbpgsmaxvec = zeros(1,maxlen);rbpgsmaxvec(1:length(rbpgsvec)) = rbpgsvec;
adimaxvec = zeros(1,maxlen);adimaxvec(1:length(adivec)) = adivec;



% display(X_jacobi);
% display(X_sdl);
% display(X_sor);
% display(X_rbpgs);
% display(X_adi);
% display(X_exact);



% Plotting the contour of the five different schemes

matrices = {X_exact, X_jacobi, X_sdl, X_sor, X_rbpgs,X_adi};
figure;
x = linspace(0, 1, n);
y = linspace(0, 1, n);
[X, Y] = meshgrid(x, y);
subplot(1, 2, 1);
contourf(X, Y, X_adi);
colorbar;
title('ADI Solution');
xlabel('X-axis');
ylabel('Y-axis');
subplot(1, 2, 2);
contourf(X, Y, X_exact);
colorbar;
title('Analytical Solution');
xlabel('X-axis');
ylabel('Y-axis');
pause(2); 


% Plotting the residual with the no. of iterations for the 5 different methods

index = log10(1:maxlen);
figure;
plot(index, log10(jacobimaxvec), 'DisplayName', 'Jacobi');
hold on;
plot(index, log10(sdlmaxvec), 'DisplayName', 'Gauss-Seidel');
plot(index, log10(sormaxvec), 'DisplayName', 'SOR');
plot(index, log10(rbpgsmaxvec), 'DisplayName', 'RBPGS');
plot(index, log10(adimaxvec), 'DisplayName', 'ADI');
xlabel('log(Iterations)');
ylabel('log(Residual)');
title('Residual vs Iterations for five different Numerical Schemes');
legend show;
hold off;


% Correct value placed again ----> Now my matrix shows the correct value of phi in x-y plane

for i = 1 : n/2
    X_jacobi([i (n-i+1)], :) = X_jacobi([(n-i+1) i], :);
    X_sdl([i (n-i+1)], :) = X_sdl([(n-i+1) i], :);
    X_sor([i (n-i+1)], :) = X_sor([(n-i+1) i], :);
    X_rbpgs([i (n-i+1)], :) = X_rbpgs([(n-i+1) i], :);
    X_adi([i (n-i+1)], :) = X_adi([(n-i+1) i], :);
    X_exact([i (n-i+1)], :) = X_exact([(n-i+1) i], :);
end    








