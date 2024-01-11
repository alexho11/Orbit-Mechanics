% Orbit Mechanics Exercise 2
% Numerical Integration of Satellite Orbit
% Hsin-Feng Ho 03770686

% Clear memory and command window
clear all; clc; close all;

% Keplerian elements of Sentinel-3
a=7191500; % semi-major axis
e=0.004; % eccentricity
I=98.3*pi/180; % inclination
Omega=257.7*pi/180; % right ascension of the ascending node
omega=144.2*pi/180; % argument of perigee
T0=0; % perigee passing time
stp=5;
% revolution period
n=sqrt(3.986004418e14/a^3);
T_rev=2*pi/n;
% calculate the cartesian coordinates of the satellite
t=0:5:3*T_rev;


%% 
[ri,ri_dot]=kep2cart(a,e,t,T0,I,Omega,omega);
y0=[ri(1,1);ri(2,1);ri(3,1);ri_dot(1,1);ri_dot(2,1);ri_dot(3,1)];
figure;
hold on;
plot3(ri(1,:),ri(2,:),ri(3,:),'LineWidth',2);
title('Orbit of Sentinel-3 in 3D','FontSize',15);
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
Earth_coast(3);
view([45 15])
saveas(gcf,'plots/3D_orbit.png');
%%
t2=0:50:3*T_rev;
[ri_50,v_50]=kep2cart(a,e,t2,T0,I,Omega,omega);

solvers={'ode45','ode23'};
option5=odeset('InitialStep',5,'MaxStep',5);
option50=odeset('InitialStep',50,'MaxStep',50);
options={option5,option50};

DR=plot_solvers(solvers, options, @yprime, t, y0,ri,ri_dot);
%%
% Plot results in radial, in-track and cross-track directions
e_r=ri./sqrt(ri(1,:).^2+ri(2,:).^2+ri(3,:).^2);
c=cross(ri,ri_dot);
e_w=cross(ri,ri_dot)./sqrt(c(1,:).^2+c(2,:).^2+c(3,:).^2);
e_s=cross(e_r,e_w);
plot_RSW(solvers,options,t,e_r,e_s,e_w,DR,@yprime);

% disturbed orbit
DR_d=plot_solvers(solvers, options, @yprime_d, t, y0,ri,ri_dot);
plot_RSW(solvers,options,t,e_r,e_s,e_w,DR_d,@yprime_d);
%%
% own integrator function
% Euler
step=5;
tspan=[0,3*T_rev];
methods = {'euler', 'rk4'};


% Loop over the methods
for j = 1:length(methods)
    method = methods{j};
    
    % Calculate the ODE
    [~, y] = odefunc(@yprime, y0, tspan, step, method);
    dr=ri-y(:,1:3)';
    dv=ri_dot-y(:,4:6)';

    % Plot the results
    figure;
    subplot(2,1,1);
    hold on;
    plot(t,dr(1,:),"LineWidth",2);
    plot(t,dr(2,:),"LineWidth",2);
    plot(t,dr(3,:),"LineWidth",2);
    ax=gca;
    ax.FontSize=13;
    ax.FontWeight="bold";
    title(['Position ',method, ' with ',num2str(step),' steps'],'FontSize',15);
    legend('x','y','z');
    xlabel('time(s)','FontSize',12,'FontWeight','bold');
    ylabel('Diffenernce(m)','FontSize',12,'FontWeight','bold');

    subplot(2,1,2);
    hold on;
    plot(t,dv(1,:),"LineWidth",2);
    plot(t,dv(2,:),"LineWidth",2);
    plot(t,dv(3,:),"LineWidth",2);
    ax=gca;
    ax.FontSize=13;
    ax.FontWeight="bold";
    title(['Velocity ',method, ' with ',num2str(step),' steps'],'FontSize',15);
    legend('x','y','z');
    xlabel('time(s)','FontSize',12,'FontWeight','bold');
    ylabel('Diffenernce(m/s)','FontSize',12,'FontWeight','bold'); 
    saveas(gcf,['./plots/',method,'_',num2str(step),'_yprime.png']);
end

step=50;
% Loop over the methods
for j = 1:length(methods)
    method = methods{j};
    
    % Calculate the ODE
    [~, y] = odefunc(@yprime, y0, tspan, step, method);
    dr=ri_50-y(:,1:3)';
    dv=v_50-y(:,4:6)';

    % Plot the results
    figure;
    subplot(2,1,1);
    hold on;
    plot(t2,dr(1,:),"LineWidth",2);
    plot(t2,dr(2,:),"LineWidth",2);
    plot(t2,dr(3,:),"LineWidth",2);
    ax=gca;
    ax.FontSize=13;
    ax.FontWeight="bold";
    title(['Position ',method, ' with ',num2str(step),' steps'],'FontSize',20);
    legend('x','y','z');
    xlabel('time(s)','FontSize',15,'FontWeight','bold');
    ylabel('Diffenernce(m)','FontSize',15,'FontWeight','bold');

    subplot(2,1,2);
    hold on;
    plot(t2,dv(1,:),"LineWidth",2);
    plot(t2,dv(2,:),"LineWidth",2);
    plot(t2,dv(3,:),"LineWidth",2);
    ax=gca;
    ax.FontSize=13;
    ax.FontWeight="bold";
    title(['Velocity ',method, ' with ',num2str(step),' steps'],'FontSize',20);
    legend('x','y','z');
    xlabel('time(s)','FontSize',15,'FontWeight','bold');
    ylabel('Diffenernce(m/s)','FontSize',15,'FontWeight','bold');
    saveas(gcf,['plots/',method,'_',num2str(step),'_yprime.png']);
end
function DR=plot_solvers(solvers, options, yprime, t, y0,r,v)
    
    % Loop over the solvers
    for i = 1:length(solvers)
        solver = str2func(solvers{i});
        
        % Loop over the options
        for j = 1:length(options)
            option = options{j};
            
            % Solve the ODE
            [t, y] = solver(yprime, t, y0, option);
            dr=r-y(:,1:3)';
            dv=v-y(:,4:6)';
            DR{i,j}=dr;

            % Check if yprime is equal to yprime_d
            if strcmp(func2str(yprime), 'yprime_d')
                title_suffix = ' disturbed';
            else
                title_suffix = '';
            end

            % Plot the results
            figure;
            subplot(2,1,1);
            hold on;
            plot(t,dr(1,:),"LineWidth",2);
            plot(t,dr(2,:),"LineWidth",2);
            plot(t,dr(3,:),"LineWidth",2);
            ax=gca;
            ax.FontSize=13;
            ax.FontWeight="bold";
            title(['Position ',func2str(solver), ' with ',num2str(option.InitialStep),' steps', title_suffix],'FontSize',20);
            legend('x','y','z');
            xlabel('time(s)','FontSize',15,'FontWeight','bold');
            ylabel('Diffenernce(m)','FontSize',15,'FontWeight','bold');

            subplot(2,1,2);
            hold on;
            plot(t,dv(1,:),"LineWidth",2);
            plot(t,dv(2,:),"LineWidth",2);
            plot(t,dv(3,:),"LineWidth",2);
            ax=gca;
            ax.FontSize=13;
            ax.FontWeight="bold";
            title(['Velocity ',func2str(solver), ' with ',num2str(option.InitialStep),' steps', title_suffix],'FontSize',20);
            legend('x','y','z');
            xlabel('time(s)','FontSize',15,'FontWeight','bold');
            ylabel('Diffenernce(m/s)','FontSize',15,'FontWeight','bold');
            saveas(gcf,['./plots/',func2str(solver),'_',num2str(option.InitialStep),'_',func2str(yprime),'.png']);
        end
    end
end


function plot_RSW(solvers,options,t,e_r,e_s,e_w,DR,yprime)
    for i = 1:length(solvers)
        solver = str2func(solvers{i});
        
        % Loop over the options
        for j = 1:length(options)
            option = options{j};
            
            d_R=dot(DR{i,j},e_r);
            d_S=dot(DR{i,j},e_s);
            d_W=dot(DR{i,j},e_w);

            % Check if yprime is equal to yprime_d
            if strcmp(func2str(yprime), 'yprime_d')
                title_suffix = ' disturbed';
            else
                title_suffix = '';
            end

            % Plot the results
            figure;
            hold on;
            plot(t,d_R,"LineWidth",2);
            plot(t,d_S,"LineWidth",2);
            plot(t,d_W,"LineWidth",2);
            ax=gca;
            ax.FontSize=13;
            ax.FontWeight="bold";
            title(['RSW ',func2str(solver), ' with ',num2str(option.InitialStep),' steps',title_suffix],'FontSize',20);
            legend('R','S','W');
            xlabel('time(s)','FontSize',15,'FontWeight','bold');
            ylabel('Diffenernce(m)','FontSize',15,'FontWeight','bold');
            saveas(gcf,['./plots/',func2str(solver),'_',num2str(option.InitialStep),'_',func2str(yprime),'_RSW.png']);
        end
    end
end