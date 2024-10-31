% Solving the advection equation:
% du/dt + c1 du/dx + c2 dudy = 0 in [-1, 1]x[-1, 1]
% BC using the linear exterapolation

clear variables; clc;

tf = 0.2;
xmin = -1; xmax = 1; ymin = -1; ymax = 1;
u1 = sqrt(2) / 2; u2 = sqrt(2) / 2;
nX = 100; 
nY = 100; % number of grid points on spatial domain x and y
x = linspace(xmin, xmax, nX); dx = x(2)-x(1);
y = linspace(ymin, ymax, nY); dy = y(2)-y(1);
dt = 0.5 / (1/dx + 1/dy); 
t = 0;

%% defining the initial phi
phi0 = zeros(nX, nY); 
xc = 0; yc = 0; R = .3;
for i=1:nX
    for j=1:nY
        phi0(i,j) = sqrt( (x(i)-xc)^2 +(y(j)-yc)^2 ) - R; 
    end
end
% we want phi equals to zero 

phin=phi0;
phinp1=zeros(nX, nY); 
while t <= tf
    if t+dt > tf
        dt = tf - t;
    end
    
%     if (u1 <= 0 && u2 <= 0)
%         for i=1:nX-1
%             for j =1:nY-1
%                 dudx = ( phin(i+1,j  ) - phin(i,j) )/dx;
%                 dudy = ( phin(i  ,j+1) - phin(i,j) )/dy;
%                 phinp1(i,j) = phin(i,j) - u1*dt*dudx - u2*dt*dudy;
%             end
%         end
%         % Treatment of boundary conditions:
%         % linear exterapolation for i = nX
%         for j=1:nY
%             phinp1(nX, j) = phinp1(nX - 1,j) + ((phinp1(nx - 1) - phinp1(nx - 2)) / dx) * (x(nX) - x(nX - 1));
%         end
%         % linear exterapolation for j = nY
%         for i=1:nX
%             phinp1(i, nY) = phinp1(i,1);
%         end
%     end
    
%     if (u1 <= 0 && u2 >  0)
%         for i=1:nX-1
%             for j =2:nY
%                 dudx = ( phin(i+1,j  ) - phin(i,j  ) )/dx;
%                 dudy = ( phin(i  ,j  ) - phin(i,j-1) )/dy;
%                 phinp1(i,j) = phin(i,j) - u1*dt*dudx - u2*dt*dudy;
%             end
%         end
%         % Treatment of boundary conditions:
%         for j=1:nY
%             phinp1(nX, j) = phinp1(1,j);
%         end
%         for i=1:nX
%             phinp1(i,  1) = phinp1(i,nY);
%         end
%     end
%     
%     if (u1 >  0 && u2 <= 0)
%         for i=2:nX
%             for j =1:nY-1
%                 dudx = ( phin(i  ,j  ) - phin(i-1,j) )/dx;
%                 dudy = ( phin(i  ,j+1) - phin(i  ,j) )/dy;
%                 phinp1(i,j) = phin(i,j) - u1*dt*dudx - u2*dt*dudy;
%             end
%         end
%         % Treatment of boundary conditions:
%         for i=1:nX
%             phinp1(i, nY) = phinp1(i,1);
%         end
%         for j=1:nY
%             phinp1(1,j) = phinp1(nX,j);
%         end
%     end
    
    if (u1 >  0 && u2 >  0)
        for i=2:nX
            for j=2:nY
                dphidx = ( phin(i  ,j  ) - phin(i-1,j  ) )/dx;
                dphidy = ( phin(i  ,j  ) - phin(i,  j-1) )/dy;
                phinp1(i,j) = phin(i,j) - u1*dt*dphidx - u2*dt*dphidy;
            end
        end
        % Treatment of boundary conditions:
        % linear exterapolation for j = 1
        for i=1:nX
            phinp1(i,1) = phinp1(i,2) + ((phinp1(i,3) - phinp1(i,2)) / dy) * (y(1) - y(2));
        end
        % linear exterapolation for i = 1
        for j=1:nY
            phinp1(1,j) = phinp1(2,j) + ((phinp1(3,j) - phinp1(2,j)) / dx) * (x(1) - x(2));
        end
    end
    
    % Update variables for the next time step:
    t = t+dt;
    phin = phinp1;
    
    ax = [min(x) max(x) min(y) max(y)];

    clf
    subplot(1,2,1), contourf(x,y,-phin,[0 0],'k-')
    axis equal, axis(ax)
    title(sprintf('geometry at t=%0.2f',t))
    subplot(1,2,2), surf(x,y,-phin,'EdgeAlpha',.2), 
    hold on
    patch([-1 1 1 -1],[-1 -1 1 1],[0 0 0 0],'k','FaceAlpha',.5)
    hold off, axp = [min(0,min(min(-phin))) max(0,max(max(-phin)))];
    axis([ax axp]), title('level set function')
    drawnow

end
% subplot(1,2,1), hold on, contourf(x,y,-phi0,[0 0],'k-')
% 
% surf(x,y,phin);
% xlim([-1 1]);
% ylim([-1 1]);
% zlim([-2 2]);
% hold on
% pause(dt*10);
% surf(x,y,phi0);