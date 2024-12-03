%*****  1D ADVECTION DIFFUSION MODEL OF HEAT TRANSPORT  *******************

%*****  Initialise Model Setup

%Image data obtained

% create x-coordinate vectors
xc = h/2:h:W-h/2;      % x-coordinate vector for cell centre positions [m]
zc = h/2:h:D-h/2;      % z-coordinate vector for cell centre positions [m]
xf = 0:h:W;            % x-coordinate vectore for cell face positions [m]
zf = 0:h:D;            % z-coordinate vectore for cell face positions [m]

[Xc,Zc] = meshgrid(xc,zc);  % create 2D coordinate arrays

%Coordinate grid on the output / for code interpretation: use this as the
%model grid 
%h   = W/Nx;
%xc  = h/2:h:W-h/2;
%zc  = h/2:h:D-h/2;
%[Xc,Zc] = meshgrid(xc,zc);

%Insulating Sides
ix3 = [           1,1:Nx,Nx      ];  
ix5 = [        1, 1,1:Nx,Nx,Nx   ];
ix7 = [   1,   1, 1,1:Nx,Nx,Nx,Nx];

%Set up Insulating Top and Bottom 
iz3 = [           1,1:Nz,Nz      ];  
iz5 = [        1, 1,1:Nz,Nz,Nz   ];
iz7 = [   1,   1, 1,1:Nz,Nz,Nz,Nz];

% create smooth random perturbation field
rng(15);
dr  = randn(Nz,Nx);

for ii = 1:10
    dr  = dr + (diff(dr(iz3,:),2,1) + diff(dr(:,ix3),2,2))/8;
end


% set initial condition for temperature at cell centres
T   = Ttop + (Tbot-Ttop)./D.*Zc + dr*1;  % initialise T array on linear gradient

% initialise density and mobility
rho = rho0.*(1 - aT.*(T-Ttop));

kT  = kT0.*ones(Nz,Nx);


% initialise output figure with initial condition
figure(1); clf
makefig(xc,zc,T,0,yr)

%*****  Solve Model Equations

dt = CFL * min((h/2),(h/2)^2/max(kT(:))); % initial time step [s]
t  = 0;  % initial time [s]
k  = 0;  % initial time step count
dTdt = 0

% loop through time steps until stopping time reached
while t <= tend

    % increment time and step count
    t = t+dt;
    k = k+1;

    % print time step header
    fprintf(1,'\n\n*****  step = %d;  dt = %1.3e;  time = %1.3e \n\n',k,dt,t)

    % store old temperature and rate
    dTdto = dTdt;
    To    = T;

    % get time step step size
    dt    = CFL * min((h/2),(h/2)^2/max(kT(:))); % time step [s]

    resnorm = 1;  % initialise residual norm
    it      = 0;  % initialise iteration count

    % loop through pseudo-transient iterations until convergence criterion reached
    while resnorm > tol

        % update temperature every 'nup' iterations
        if ~mod(it,nup) && k>=1
            ssssss = 31546432
            disp(ssssss)
            % get T-dependent segregation mobility
            kT  = kT0 + cT.*max(0,T).^mT;

            % get rate of change
            dTdt = diffusion(T,kT,h,ix3,iz3);

            % get temperature residual
            res_T = (T - To)/dt - (dTdt + dTdto)/2;

            % set isothermal boundaries on top/bot
            res_T(1  ,:) = 0;
            res_T(end,:) = 0;

            % get solution update
            upd_T = - alpha*res_T*dt/2;

            % update solution
            T     = T + upd_T;

        end

        % get density and density contrast
        rho   = rho0.*(1 - aT.*(T-Ttop));  % T-dependent density
        Drho  = rho - mean(rho,2);         % subtract horizontal mean
        Drhoz = (Drho(iz3(1:end-1),:)+Drho(iz3(2:end),:));  % on z-faces
        Drhoz([1 end],:) = 0;              % no flow across top/bot bounds

        it = it+1; % increment iteration count

        % get residual norm and print convergence every 'nup' iterations
        if ~mod(it,nup)
            resnorm = norm(upd_T(:),2)./norm(T(:)+eps,2) ...
                    + norm(upd_p(:),2)./norm(p(:)+eps,2);
            if isnan(resnorm); error('!!! Solver failed with nan !!!'); end
            fprintf(1,'     it = %d;  res = %e \n',it,resnorm); 
        end

    end

    % plot model progress every 'nop' time steps
    if ~mod(k,nop)
        makefig(xc,zc,T,t,yr);
    end

end


%*****  Utility Functions  ************************************************

% Function to make output figure
function makefig(x,z,T,t,yr)

clf; 

% plot temperature in subplot 1
subplot(2,2,1);
imagesc(x,z,T); axis equal tight; colorbar; hold on
contour(x,z,T,[100,150,200],'k');

ylabel('z [m]','FontSize',15)
title('Temperature [C]','FontSize',17)

end

% Function to calculate diffusion rate
function [dTdt] = diffusion(f,k,h,ix,iz)
% calculate heat flux by diffusion

qx = - k .* diff(f(:,ix), 1, 2)/h;
qz = - k .* diff(f(iz, :), 1, 1)/h;

% calculate flux balance for rate of change
dTdt = - (diff(qx)/h+diff(qz)/h);


end


