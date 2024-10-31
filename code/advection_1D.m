% Solving the advection equation in 1D:
% du/dt + c du/dx = 0 in [-1, 1]
% uInitial = piecewise function
% BC using the linear exterapolation

tf = 0;
xmin = -1; xmax = 1;
u1 = 1;
u2 = -1;
nGridPoints = 100;
x = linspace(xmin, xmax, nGridPoints); dx = x(2)-x(1);
unp1 = zeros(nGridPoints, 1);
dt = 0.5*dx;
t = 0;

% initial function 
xl = -1/3; % denote as x-low 
xm = 1/3; % denote as x-max
u0 = zeros(nGridPoints, 1);
for i = 1:nGridPoints 
    if i <= floor(((xl - xmin)/dx ) + 1) 
        u0(i) = 0;
    elseif floor(((xl - xmin)/dx ) + 1) <= i && i <= floor(((xm - xmin)/dx ) + 1)
        u0(i) = 1; 
    else
        u0(i) = 0;
    end 
end

un = u0;
un2 = u0; 
unp12 = zeros(nGridPoints, 1); 

while t <= tf
    if t+dt > tf
        dt = tf - t;
    end
    
    if (u1 <= 0)
        % Update the inner grid points:
        for i=1:nGridPoints-1
            dudx = ( un(i+1) - un(i) )/dx;
            unp1(i) = un(i) - u1*dt*dudx;
        end
        % Impose the periodic boundary condition:
        unp1(nGridPoints) = unp1(nGridPoints - 2) + ((unp1(nGridPoints - 1) - unp1(nGridPoints - 2)) / dx) * (x(nGridPoints) - x(nGridPoints - 2));
    else % u1 > 0
        for i=2:nGridPoints
            dudx = ( un(i) - un(i-1) )/dx;
            unp1(i) = un(i) - u1*dt*dudx;
        end
        % Impose the periodic boundary condition:
        unp1(1) = unp1(2) + ((unp1(3) - unp1(2)) / dx) * (x(1) - x(2));
    end
    
    % Update variables for the next time step:
    t = t+dt;
    un = unp1;
    
    if (u2 <= 0)
        % Update the inner grid points:
        for i=1:nGridPoints-1
            dudx = ( un2(i+1) - un2(i) )/dx;
            unp12(i) = un2(i) - u2*dt*dudx;
        end
        % Impose the periodic boundary condition:
        unp12(nGridPoints) = unp12(nGridPoints - 2) + ((unp12(nGridPoints - 1) - unp12(nGridPoints - 2)) / dx) * (x(nGridPoints) - x(nGridPoints - 2));
    else % u2 > 0
        for i=2:nGridPoints
            dudx = ( un2(i) - un2(i-1) )/dx;
            unp12(i) = un2(i) - u2*dt*dudx;
        end
        % Impose the periodic boundary condition:
        unp12(1) = unp12(2) + ((unp12(3) - unp12(2)) / dx) * (x(1) - x(2));
    end
    
    % Update variables for the next time step:
    t = t+dt;
    
    un2 = unp12;
    plot(x, u0, 'b:', x, unp12, 'g-.', x, unp1, '--r', 'LineWidth',2.25);  legend('initial phi', 'negative velocity', 'positive velocity'); hold on; 
    s = sprintf('phi at t=%2.2f', t)
    title(s); xlabel('x'); ylabel('phi values'); 
%     pause(dt); hold off; 
%     fixfig; 

    if t == tf
        break
    end
    
end

