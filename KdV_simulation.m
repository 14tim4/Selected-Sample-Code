% This numerically solves the KdV equation using fft and for time-stepping, using Crank-Nicolson. See the readme file for further description.
% Note some minor changes have been made from the original code used in the assigment in 2020.
%
% Author: Timothy Lapuz

clear;
clc;



% Initialisations
n = 1024;
l = 100;
dx = l/(n-1);
x = 0:dx:l;
dk = 2*pi/(n*dx); % Fourier wave numbers
%k = [0:dk:n/2*dk,-(n/2-1)*dk:dk:-dk]; 
k = [-(n/2-1)*dk:dk:n/2*dk]; % Fourier wave numbers array
%k = [0:dk:(n-1)*dk];


% Time
dt = 0.0025;
tmax = 100; 
nmax=round(tmax/dt);
nplots = 200;
nsnap=round(nmax/nplots);


% Linear and nonlinear operators
L = 1i*k.^3;
N = -3*1i*k;


% We set the paramaters for the ICs
cover = l/8; % We want to cover 1/8 of the domain with the Gaussian curve
sigma = cover/6; % We aim for 1/8 of the domain to be within 3 standard deviations
A = [0.5 1 1.5 2 2.5];
B = 1/(2*sigma^2); % Calculate B


for p = 1:length(A)
    u0 = A(p)*exp(-B*(x - l/2).^2);

    % Record the first one
    uplot(1,:) = u0;

    % Fourier transform of IC
    U0 = fftshift(fft(u0));
    U0sq = fftshift(fft(u0.^2));

    % Initialise array for time
    t(1)=0;

    % Euler for the first timestep
    k1 = dt*(L.*U0 + N.*U0sq);
    U1 = U0 + k1; % Euler integration in k-space

    % Inverse Fourier
    u1 = real(ifft(ifftshift(U1)));

    % Initialise array for time
    t(2)=dt;

    % Crank-Nicolson Method for the rest of the time-stepping
    for m=1:nplots
        for j=1:nsnap

            % Need the Fourier transforms
            U0sq = fftshift(fft(u0.^2));
            U1sq = fftshift(fft(u1.^2));
            U0 = fftshift(fft(u0));
            U1 = fftshift(fft(u1));

            % Time integration in k-space
            U2 = (1+dt/2*L).*U1./(1-dt/2*L) + dt*(1.5*N.*U1sq-0.5*N.*U0sq)./(1-dt/2*L);

            % Inverse Fourier
            u0 = real(ifft(ifftshift(U1)));
            u1 = real(ifft(ifftshift(U2)));
        end

      % We record the current terms
      U = U2;

      % Inverse Fourier
      u = real(ifft(ifftshift(U)));

      % Record the values for a snapshot and update the time.
      t(m+1)=t(m)+nsnap*dt;
      uplot(m+1,:)=u;

      % Update waitbar and message
      waitbar(m/nplots);
    end

    if p == 1
        uplot1 = uplot;
    elseif p == 2
        uplot2 = uplot;
    elseif p == 3
        uplot3 = uplot;
    elseif p == 4
        uplot4 = uplot;
    else
        uplot5 = uplot;
    end
        

end
      
%% Plotting
uplot11 = uplot1(31,:);
uplot22 = uplot2(31,:);
uplot33 = uplot3(31,:);
uplot44 = uplot4(31,:);
uplot55 = uplot5(31,:);

% Some plots
figure(1)
plot(x,uplot11,'Color',[0 0 1],'LineWidth', 1);
hold on
plot(x,uplot22,'Color',[0 1 0],'LineWidth', 1);
hold on
plot(x,uplot33,'Color',[1 0 0],'LineWidth', 1);
hold on
xlabel('x','FontSize', 12);
ylabel('u', 'FontSize', 12);
title('Plot at $t=15$','Interpreter','latex', 'FontSize', 16)
ylim([0 4]);
legend('A = 0.5','A = 1','A = 1.5');
set(gca,'FontSize',16);


figure(2)
plot(x,uplot44,'Color',[0 0 1],'LineWidth', 1);
hold on
plot(x,uplot55,'Color',[0 1 0],'LineWidth', 1);
hold on
ylim([0 4]);
xlabel('x','FontSize', 12);
ylabel('u', 'FontSize', 12);
title('Plot at $t=15$','Interpreter','latex', 'FontSize', 16)
legend('A = 2','A = 2.5');
set(gca,'FontSize',16);

