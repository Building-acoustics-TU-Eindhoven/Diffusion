%% Diffusion equation for sound field prediction
% 2D implementation
% Based on Navarro Applied Acoustics (2012)

clc
clear all

roomname1 = 'QR_LR_02'; name = 'Room';
filename4 = [roomname1,'_results_01.mat'];
% filename4 = [roomname1,'_results_street_01.mat'];
filename4 = [roomname1,'_results_street_02.mat'];
load (filename4,'Dx_3oct_rcv_meanX_fitted')
Dy_wbm1 = 55.5;
Dy_wbm2 = 46.46;

c = 343;
rho=1.21;

% th type of constant A on boundary conditions . options:
th = 1;
if th==1
    Acoef = @(alpha) (c*alpha)/4; % Sabine
elseif th==2
    Acoef = @(alpha) (c*(-log(1-alpha)))/4; % Eyring
elseif th==3
    Acoef = @(alpha) (c* alpha)/(2*(2-alpha)); % Xiang
end

alphax = 0.99;
alphay = 0.1;
Absx = Acoef(alphax);
Absy = Acoef(alphay);

dx = 0.1;    % Spatial discretization [m]

dt = 0.0001;    % Temporal discretization (time step) [s]
fm = 1/dt;

% Room dimensions
lx = 39.7; % Length [m]
ly = 2.8; % Width [m]

S2D = lx*ly;  % Surface [m2]
P2D = 2*(lx+ly); % Perimeter [m]

x = 0:dx:lx;
y = 0:dx:ly;

Nx = length(x);
Ny = length(y);

[yy,xx] = meshgrid(y,x)


% Diffusion parameters
m_atm = 0;   % Atmospheric attenuation
mfp = pi*S2D/P2D;    % Mean-free-path
%Dx = ones(size(xx)).*mfp*c/2;
%Dy = ones(size(xx)).*mfp*c/2;
Dx = mfp*c/3;    % Diffusion coefficient
%Dx = ones(size(xx)).*mfp*c/2;
% Dy = ones(size(xx)).*Dy_wbm1;
%Dy = ones(size(xx)).*Dy_wbm2;

%Dx_vec = abs([Dx_3oct_rcv_meanX_fitted Dx_3oct_rcv_meanX_fitted(end) Dx_3oct_rcv_meanX_fitted(end) Dx_3oct_rcv_meanX_fitted(end) Dx_3oct_rcv_meanX_fitted(end)]');
%Dx = repmat(Dx_vec,1,Ny);
%inxDx = find(Dx<0);


% Beta condition
beta_zerox = (2*Dx*dt)/(dx^2);
beta_zeroy = (2*Dx*dt)/(dx^2);
beta_zero = beta_zerox + beta_zeroy;

stabcond = (beta_zero.^2-1)./(1 + beta_zero.^2 + 2.*beta_zero);

if max(stabcond(:)) > 1
    return
end

% Initial conditions (excitation with Gaussian)
%source density energy  flux 
Ws=10^-2;   % Point power
Vs=4/3*pi*dx^3; % Source volume
w1 = Ws/Vs;

recording_time = 3;
recording_steps = ceil(recording_time/dt);
sourceon_time = 1;% Interrupted source
sourceon_steps = ceil(sourceon_time/dt);

s1= w1.*ones(1,sourceon_steps);               
source = [s1 zeros(1,recording_steps-sourceon_steps)];      %for source on time. interrupted source

x_sou = 0.5;
y_sou = 0.8;
indx_sou = dsearchn([xx(:), yy(:)],[x_sou y_sou]);

x_rec = [10];
y_rec = [1.5];
for ii = 1:length(x_rec)
    indx_rec(ii) = dsearchn([xx(:), yy(:)],[x_rec(ii) y_rec(ii)]);
end


s = single(zeros(Nx,Ny));
s(indx_sou) = source(1);


w_new = zeros(Nx,Ny);
w = w_new;
w_old = w;


for steps = 1:recording_steps;
    time = steps*dt
    s(indx_sou) = source(steps);

    wiminus1 = w([1 1:Nx-1],[1:Ny]);
    wiplus1 = w([2:Nx Nx],[1:Ny]);
    wjminus1 = w([1:Nx],[1 1:Ny-1]);
    wjplus1 = w([1:Nx],[2:Ny Ny]);
    
    w_new = (w_old.*(1-beta_zero) - 2*dt*c*m_atm*w + 2*dt*s +...
        beta_zerox.*(w([2:Nx Nx],[1:Ny]) + w([1 1:Nx-1],[1:Ny])) +...
        beta_zeroy.*(w([1:Nx],[2:Ny Ny]) + w([1:Nx],[1 1:Ny-1])))./(1+beta_zero);
    
    w_new(1,:) = (4*w_new(2,:)-w_new(3,:))./(3+(2*Absx*dx./Dx(1,:)));
    w_new(end,:) = (4*w_new(end-1,:)-w_new(end-2,:))./(3+(2*Absx*dx./Dx(end,:)));
    w_new(:,1) = (4*w_new(:,2)-w_new(:,3))./(3+(2*Absy*dx./Dx(:,1)));
    w_new(:,end) = (4*w_new(:,end-1)-w_new(:,end-2))./(3+(2*Absy*dx./Dx(:,end)));
    
    w_old = w;
    w = w_new;
    
    w_rec(steps,:) = w_new(indx_rec);
       
%     figure(1)
%     surf(xx,yy,abs(w_new))    
%     colormap hot
%     view(2)
%     light;
%     lighting phong;
%     camlight('left');
%     axis equal
%     shading interp;
%     drawnow

end
    
figure(2)
plot((1:steps)*dt,10*log10(abs(w_rec)*rho*c^2/(2e-5)^2))
hold on

figure(3)
plot(xx(:,14),10*log10(abs(w_new(:,14))*rho*c^2/(2e-5)^2))
hold on


