 %************   GEOTHERMAL PROJECT (VERS. B)  *************



%*** MODEL SETUP *************************************************

% creating grid of cells and setting position for centre or face
x_cc = dx/2:dx:W-dx/2; % x coord: cell centres to start at centre of first cell of width dx and continue in matrix along 
x_fc = 0:dx:W; % x cord: cell faces set as matrix entries with boundaries at each edge of cell.
z_cc = dx/2:dx:D-dx/2; % z coord: same thing but depth not width 
z_fc = 0:dx:D; % z coord: same thing but depth not width

[Xc,Zc] = meshgrid(x_cc,z_cc);  % use the 1D matrices of cell centres to create 2D matrix of domain

% BOUNDARY CONDITIONS: add index cells on each side to apply
ix = [ 1,1:Nx,Nx ];  % closed/insulating sides
iz = [ 1,1:Nz,Nz ];  % closed/insulating top, flux grad at bottom


% INITIAL CONDITIONS: temperature
T   = T0 + dTdz_0(2).*Zc;  % initialise T array on linear gradient
Ta = T; % store this initial temperature
air = units == 4; % name air if the unit is 4 (depending on image)

u = zeros(Nz, Nx);
w = zeros(Nz, Nx);

%*** MODEL EQUATIONS *************************************************

t = 0;  % initial time [s]
tau = 0;  % initial time step count
dt = CFL * min((dx/2)^2 / max(k0(:)), dx / max(abs(w(:)) + eps)); % finding the time-step using the max ko value, spatial step and CFL [s]

while t <= tend % for each t until the end time has been reached

    % increment time and step count
    t = t+dt; % counting time
    tau = tau+1; % counting step

    % reset air section to be fixed value
    T(air) = Tair;
    
    % create a temperature depdent density for the pressure gradient 
    rho_P = rho0.* (1 - aT.* (T-T0));

    % calculate the hydrostatic pressure (pressure gradient)
    P = zeros(Nz,Nx);
    for i = 2:Nz 
        P(i,:) = P(i-1,:) - rho_P(i,:) * g * dx;
    end
    
    % Darcy velocity from P gradient
    w = -(k_p ./ mu) .* [diff(P,1,1); zeros(1,Nx)]/ dx;

    w = [w(1,:);(w(1:end-1,:) + w(2:end,:))/2];


            % 4th-order Runge-Kutta time integration scheme
                    
            dTdt1 = diffusion(T,dTdz_0,            k0,dx,ix,iz) - advection(T,u,w,dx,ix,iz);
            dTdt2 = diffusion(T+dTdt1/2*dt, dTdz_0,k0,dx,ix,iz) - advection(T+dTdt1/2*dt,u,w,dx,ix,iz);
            dTdt3 = diffusion(T+dTdt2/2*dt, dTdz_0,k0,dx,ix,iz) - advection(T+dTdt2/2*dt,u,w,dx,ix,iz);
            dTdt4 = diffusion(T+dTdt3  *dt, dTdz_0,k0,dx,ix,iz) - advection(T+dTdt3*dt,u,w,dx,ix,iz);
        
            T = T + (dTdt1 + 2*dTdt2 + 2*dTdt3 + dTdt4)/6 * dt + (Hr ./ rho ./ Cp)*dt;
  
   
    if ~mod(tau,nop) || t >= tend
    makefig(x_cc, z_cc, T, t, yr); % call the plotting function

         if t >= tend
             % Save figure at final time step
             filename = ['Temperature_Final_', num2str(round(t/yr)), '_yrs.png'];
             saveas(gcf, filename);  % saves as PNG
         end
    end
end



%*** MODEL FUNCTIONS *************************************************

function [dTdt] = diffusion(f,dTdz_0,k0,dx,ix,iz) % diffusion rate function

% average k0 values to get cell face values
kx = (k0(:, ix(1:end-1)) + k0(:, ix(2:end)))/2; % averaged diffusion coefficient kx in the x-direction
kz = (k0(iz(1:end-1), :) + k0(iz(2:end), :))/2; % averaged diffusion coefficient kz in the z-direction

% calculate heat flux by diffusion
qx = - kx .* diff(f(:, ix), 1, 2)/dx;
qz = - kz .* diff(f(iz, :), 1, 1)/dx;

% set boundary conditions - using geothermal gradient at bottom boundary
qz(end, :) = -kz(end, :) .* dTdz_0(2);

% calculate flux balance for rate of change
dTdt = -(diff(qx, 1, 2)/dx + diff(qz, 1, 1)/dx);

end

function [dTdt] = advection(T,u,w,dx,ix,iz)

Tx = (T(:, ix(1:end-1)) + T(:, ix(2:end))) / 2;
Tz = (T(iz(1:end-1), :) + T(iz(2:end), :)) / 2;

ux = (u(:, ix(1:end-1)) + u(:, ix(2:end))) /2;
wz = (w(iz(1:end-1),:) + w(iz(2:end),:)) /2;

advx = -diff(ux.* Tx,1,2)/dx; %divergence of fluxes
advz = -diff(wz.*Tz,1,1)/dx;

dTdt = advx + advz;

end

% Function to make output figure
function makefig(x,z,T,t,yr)

clf; 

%plot temperature across domain
imagesc(x,z,T); axis equal tight; colorbar; hold on
ylabel('z [m]','FontSize',18,'FontName','Times New Roman')
xlabel('x [m]','FontSize',18,'FontName','Times New Roman')
title(['Temperature; time = ',num2str(t/yr),'yr'],'FontSize',20,'FontName','Times New Roman');
ylabel(colorbar, 'Temperature [°C]','FontSize',18,'FontName','Times New Roman');  
[C,dx] = contour(x,z,T, [40,90, 150],'r','Linewidth',2);
clabel(C,dx,'Fontsize',15,'Color','r','FontName','Times New Roman');

drawnow;

end 
