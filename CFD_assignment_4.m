clc;clear all;
%% Initialization
L = 1;                                      %Length
W = 1;                                      %Height or width
dx = 0.1;
Nx = floor(L/dx);
Ny = floor(W/dx);

T_analytical(1:Nx,1:Ny) = 0;                %Initialization of variables
T1 = 100;
T2 = 200;
n = 100;

%% Analytical
for i = 0:Ny                                %for y axis/row of Temp matrix
    for j = 0:Nx                            %for x axis/column of Temp matrix
        
        x = j*dx;
        y = i*dx;

       s = 0;
       for k = 1:n
             s = s + (((((-1)^(k+1))+1)/k)*sin(k*pi*x/L)*sinh(k*pi*y/L)/sinh(k*pi*W/L)); %Analytical Formula
       end
       G = (2/pi)*s;
       T_analytical(i+1,j+1) = T1 + (T2 - T1)*G;
    end
end


%% Numerical
Nx = 10;                %Put the Number of nodes here (x direction)
Ny = 10;                %Put the Number of nodes here (y direction)
T(1:Nx+2,1:Ny+2) = 150;
dx = 1/Nx;
dy = 1/Ny;      

%Alpha Calculation
rho = 7750;
cp = 500;
k = 16.2;
alpha = k/(rho*cp);

%Time Step Calculation
c = 0.4;
dt = (c*dx^2)/(alpha*(1+(dx/dy)^2));

No_of_iteration = 1000;

for iter = 1:No_of_iteration
    
    T_old = T;           %Temp at nth iteration
    
    %Updating Boundary Condition (fictious cells)
    T(1,:) = 2*T1-T(2,:);                  %Bottom Layer
    T(Ny+2,:) = 2*T2-T(Ny+1,:);            %Top Layer  
    T(:,1) = 2*T1-T(:,2);                  %Right layer
    T(:,Nx+2) = 2*T1-T(:,Nx+1);            %Left layer
    T_new = T;

    
    %Internal Cells
    for i = 2:Nx+1
        for j = 2:Ny+1
            T_new(i,j) = (1-(2*alpha*dt/dx^2)-(2*alpha*dt/dy^2))*T(i,j) + (alpha*dt/dx^2)*(T(i-1,j) + T(i+1,j)) + (alpha*dt/dy^2)*(T(i,j+1) + T(i,j-1));
        end
    end
    T = T_new;                     %Temp at (n+1)th iteration
    
    %%Taking the Central Cell
    Central_temp(iter) = T(int8(Nx/2)+1,int8(Ny/2)+1);
    
    %%Max temp in matrix for grid covergence
    temp = T(2:Nx+1,2:Ny+1);
    maxtemp = max(temp(:));
    
    %%Error estimation
    err = T_new - T_old;
    err = err(2:Nx+1,2:Ny+1);
    err = err.*err;
    err = sqrt(sum(err(:))/(Nx*Ny));
    %err = abs(sum(err(:))/(Nx*Ny));
    if err < 1e-6 && err ~= 0
        fprintf('Steady state Convergence is reached .....!!!!')
        break
    end
    
end
T_numerical = T;

%% Plotting
figure(1)
contourf(T_analytical,20)                        %Plotting
title('Analytical Temperature')
colorbar

figure(2)
contourf(T_numerical(2:Nx+1,2:Ny+1),20)
title('Numerical Temperature')
colorbar

%%Center point temperature vs time
t = 1:iter;
figure(3)
plot(t,Central_temp)
title('Center point temperature vs time')
%% Taking Temp Data from x = 0.5 and y = 0.5
x = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
node_x = int8(x/dx);
node_y = int8(x/dy);
temp = T_numerical(2:Nx+1,2:Ny+1);
for i = 1:length(node_y)
    yline(i) = temp(node_x(i),int8(0.5/dx));
end
for i = 1:length(node_x)
    xline(i) = temp(int8(0.5/dx),node_y(i));
end


%% Heat Flux Calulation(Global Conservation)
q = -k*(T_numerical(2,:)-T_numerical(1,:))/dy;     %Bottom layer
q = q(2:Nx+1);
q_bottom = sum(q)

q = -k*(T_numerical(Ny+1,:)-T_numerical(Ny+2,:))/dy;     %Top layer
q = q(2:Nx+1);
q_top = sum(q)

q = -k*(T_numerical(:,Nx+1)-T_numerical(:,Nx+2))/dx;     %Right layer
q = q(2:Ny+1);
q_right = sum(q)   

q = -k*(T_numerical(:,2)-T_numerical(:,1))/dx;     %Left layer
q = q(2:Ny+1);
q_left = sum(q) 

residue = q_left + q_right + q_top + q_bottom