function video_fun(x,thetaM,torque,vrel,titledescription,params,record,framerate,videoname)
%% curves_video_fun
% Levi Manring, Duke University
% 2021
%
% This function is used to provide a video demonstration of the vehicle
% moving on the surface profile. The upper plot shows a visualization of
% the vehicle moving on the surface and the lower plot shows a
% visualization of the non-dimensional applied torque to wheel A and
% the non-dimensional relative velocity of wheel A. Inputs torque and vrel
% can be changed for any other useful metrics the user is interested in
% visualizing during the course of a simulation.
%
% Inputs:
%   x: Nx1 array containing state x
%   torque: Nx1 array containing the applied torque to wheel A
%   vrel: Nx1 array containing the relative velocity of wheel A with
%           respect to the surface profile
%   titledescription: a string containing the desired title for the plot
%   params: a parameter structure including the masses, moments of
%           inertia, and dimensions needed to define the planar vehicle model
%   curves: a curves structure including the interpolated dynamics model
%           parameters needed to define the planar vehicle model
%   record: 1x1 boulean which is used to indicate whether or not the user
%           wishes to record the video file (1 means record, 0 means do not record)
%   framerate: 1x1 double indicating the framerate to record the video file
%   videoname: a string containing the desired video file name

%% Initialize

% Create the figure and set it to be the screen size
movie_fig = figure; hold on; 
set(gcf, 'Position', get(0, 'Screensize'));
clf(1);

% Unpack the parameters structure
l = params.dim.l;
R = params.dim.R;
h = params.dim.h;
xc = params.dim.xc;
yc = params.dim.yc;
dt = params.dt;

% Rescale vrel and torque so they fit on the same plot
if abs(max(vrel)) ~= 0
    vrel_sc = vrel/max(abs(vrel));
else
    vrel_sc = vrel;
end
torque_sc = torque/max(abs(torque));

% Initialize plotting vectors
th = 0:pi/50:2*pi; % radians vector for wheels
Lx = length(x);
Y = zeros(Lx,1);
x_wheel_A    = zeros(Lx,length(th));
y_wheel_A    = zeros(Lx,length(th));
x_wheel_B    = zeros(Lx,length(th));
y_wheel_B    = zeros(Lx,length(th));
x_mass_M     = zeros(Lx,5);
y_mass_M     = zeros(Lx,5);
x_cm_M_lines = zeros(Lx,3);
y_cm_M_lines = zeros(Lx,3);
x_cm_M_point = zeros(Lx,1);
y_cm_M_point = zeros(Lx,1);

% Generate the plotting information for the vehicle shape given state x
for k = 1:Lx
    xk = x(k,1);
    
    % Get surface profile, its derivative and vehicle angle
    [Y(k,1), Yx, ~, ~, ~] = y_eval_fcn(xk,params);
    theta_m = thetaM(k,1);
    
    % Get position vectors for wheel A and wheel B
    R_A_i = xk - R*Yx/sqrt(1 + Yx^2);
    R_A_j = Y(k,1) + R/sqrt(1 + Yx^2);
    R_B_i = R_A_i + l*cos(theta_m);
    R_B_j = R_A_j + l*sin(theta_m);
    
    % create x y coordinates to plot wheel A and B
    x_wheel_A(k,:) = R*cos(th) + R_A_i;
    y_wheel_A(k,:) = R*sin(th) + R_A_j;
    x_wheel_B(k,:) = R*cos(th) + R_B_i;
    y_wheel_B(k,:) = R*sin(th) + R_B_j;
    
    % create x y coordinates to plot body mass M (lc = left corner, rc = right corner)
    xlc = R_A_i - h*sin(theta_m);
    ylc = R_A_j + h*cos(theta_m);
    xrc = xlc + l*cos(theta_m);
    yrc = ylc + l*sin(theta_m);
    x_mass_M(k,:) = [R_A_i xlc xrc R_B_i R_A_i];
    y_mass_M(k,:) = [R_A_j ylc yrc R_B_j R_A_j];
    
    % create center of mass lines for body mass M (cm = center of mass)
    x_cm1 = R_A_i - yc*sin(theta_m);
    y_cm1 = R_A_j + yc*cos(theta_m);
    x_cm2 = x_cm1 + xc*cos(theta_m);
    y_cm2 = y_cm1 + xc*sin(theta_m);
    x_cm3 = R_A_i + xc*cos(theta_m);
    y_cm3 = R_A_j + xc*sin(theta_m);
    x_cm_M_lines(k,:) = [x_cm1 x_cm2 x_cm3];
    y_cm_M_lines(k,:) = [y_cm1 y_cm2 y_cm3];
    x_cm_M_point(k,1) = x_cm2;
    y_cm_M_point(k,1) = y_cm2;
    
end

% Generate refined plot of the surface profile
X_refine = linspace(min(x)-params.dim.l,max(x)+params.dim.l,12*length(x))';
Y_refine = zeros(size(X_refine));
for m = 1:length(X_refine)
    [Y_refine(m,1), ~, ~, ~, ~] = y_eval_fcn(X_refine(m,1),params);
end

%% Create Upper Plot
upper = subplot(2,1,1);
upper.XTick = []; upper.YTick = []; 
upper.YAxis.Visible = 'off'; upper.XAxis.Visible = 'off';
upper.TitleFontSizeMultiplier = 2;
title(titledescription,'Interpreter','Latex');
grid off; hold on;
p_surfaceprofile = plot(X_refine,Y_refine,'m-','linewidth',2);
p_wheel_A       = plot(x_wheel_A(1,:),y_wheel_A(1,:),'k-','linewidth',2);
p_wheel_B       = plot(x_wheel_B(1,:),y_wheel_B(1,:),'k-','linewidth',2);
p_mass_M        = plot(x_mass_M(1,:),y_mass_M(1,:),'b-','linewidth',2);
p_cm_M_lines    = plot(x_cm_M_lines(1,:),y_cm_M_lines(1,:),'g-','linewidth',2);
p_cm_M_point    = plot(x_cm_M_point(1,:),y_cm_M_point(1,:),'ro','MarkerFaceColor','r');

set(upper,'fontsize',18,'TickLabelInterpreter','Latex');
upper.XLim = [x(1,1)-32*R, x(1,1)+32*R];
upper.YLim = [Y(1,1)-12*R, Y(1,1)+12*R];

%% Create Lower Plot
lower = subplot(2,1,2);
grid on; hold on;
p_vrel        = plot(0,vrel_sc(1,1),'b-','linewidth',2);
p_torque      = plot(0,torque_sc(1,1),'r-','linewidth',2);
p_vrel_inst   = plot(0,vrel_sc(1,1),'bo','markersize',16,'linewidth',2);
p_torque_inst = plot(0,torque_sc(1,1),'rx','markersize',16,'linewidth',2);
xlabel('$t$ (s)','Interpreter','Latex'); ylabel('Nondimensional Values','Interpreter','Latex');
legend('$\hat{v}_{rel}$','$\hat{\tau}_A$','Interpreter','Latex'); 
set(lower,'fontsize',18,'TickLabelInterpreter','Latex');
lower.XLim = [0, dt*(k-1)];
lower.YLim = [-1, 1];

frames = struct('cdata',[],'colormap',[]);
frames(1) = getframe(gcf);

%% Update the plot data for each frame
for k = 2:length(x)
    p_wheel_A.XData = x_wheel_A(k,:);
    p_wheel_A.YData = y_wheel_A(k,:);
    
    p_wheel_B.XData = x_wheel_B(k,:);
    p_wheel_B.YData = y_wheel_B(k,:);
    
    p_mass_M.XData = x_mass_M(k,:);
    p_mass_M.YData = y_mass_M(k,:);
    
    p_cm_M_lines.XData = x_cm_M_lines(k,:);
    p_cm_M_lines.YData = y_cm_M_lines(k,:);
    
    p_cm_M_point.XData = x_cm_M_point(k,:);
    p_cm_M_point.YData = y_cm_M_point(k,:);
    
    p_vrel.XData = [p_vrel.XData, dt*(k-1)];
    p_vrel.YData = [p_vrel.YData, vrel_sc(k,1)];
    
    p_torque.XData = [p_torque.XData, dt*(k-1)];
    p_torque.YData = [p_torque.YData, torque_sc(k,1)];
    
    p_vrel_inst.XData = dt*(k-1);
    p_vrel_inst.YData = vrel_sc(k,1);
    
    p_torque_inst.XData = dt*(k-1);
    p_torque_inst.YData = torque_sc(k,1);
    
    upper.XLim = [x(k,1)-32*R, x(k,1)+32*R];
    upper.YLim = [Y(k,1)-12*R, Y(k,1)+12*R];
    
%     pause(params.dt);
    % update the plot now and save the frame
    drawnow
    frames(k) = getframe(gcf);
    
end

% Choose whether or not to record the video to a file
if record == 1
    
    v = VideoWriter(videoname,'MPEG-4');
    v.FrameRate = framerate;
    open(v);
    writeVideo(v,frames);
    close(v);
    
end






