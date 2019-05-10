function [vel, omg] = vel_estimation(sensor, varargin)

persistent pointTracker;
persistent oldXc;
persistent oldYc;
persistent oldZc;
persistent oldSensor;
persistent invK;
persistent cam_invK;
persistent framei;
persistent olddt;
persistent oldv;
persistent oldomg;

N_tracking = 150; 

if isempty(oldSensor) | oldSensor.t > sensor.t
    % initialize last frame
    framei = 1;
    olddt = 0.02;
    oldv = zeros(3,1);
    oldomg = zeros(3,1);
    % init KLT
    pointTracker = vision.PointTracker('MaxBidirectionalError', 2);
    corners = detectFASTFeatures(sensor.img);
    points = corners.selectStrongest(N_tracking);
    points = double(points.Location);
    initialize(pointTracker, points, sensor.img);
    
    % init constants
    K = [311.0520 0 201.8724; 0 311.3885 113.6210; 0 0 1];
    invK = inv(K);
    cam_K = [311.0520 0 201.8724; 0 311.3885 113.6210; 0 0 1];
    cam_invK = inv(cam_K);
    
    % estimate poses
    [Rc_w, Tc_w] = estimate_cam_pose(sensor, cam_invK);
    [Xc, Yc, Zc, ~] = get_XYZc(Rc_w, Tc_w, points, invK);
    
    % backup
    oldXc = Xc;
    oldYc = Yc;
    oldZc = Zc;
    oldSensor = sensor;
    
    vel = zeros(3,1);
    omg = zeros(3,1);
    return
end

% tic
points = step(pointTracker, sensor.img);
% toc
[h,w]=size(sensor.img);
ids = points(:,1) > w*0.12 & points(:,1) < w*0.88 & points(:,2) > h*0.12 & points(:,2) < h*0.88;
points = points(ids,:);

% estimate pose
[Rc_w, Tc_w] = estimate_cam_pose(sensor, cam_invK);
[Xc, Yc, Zc, ~] = get_XYZc(Rc_w, Tc_w, points, invK);
alpha = 0.05;
dt = sensor.t - oldSensor.t;
x = Xc./Zc;
y = Yc./Zc;
oldx = oldXc(ids)./oldZc(ids);
oldy = oldYc(ids)./oldZc(ids);

% solve for the velocities
dot_x = [x - oldx]'./dt;
dot_y = [y - oldy]'./dt;
A_all = get_A(x',y',Zc');
v_omega_c = A_all \ [dot_x;dot_y];

% redetect the points
if mod(framei, 2) == 0
    corners = detectFASTFeatures(sensor.img);
    points = corners.selectStrongest(N_tracking);
    points = double(points.Location);

    [h,w]=size(sensor.img);
    ids = points(:,1) > w*0.12 & points(:,1) < w*0.88 & points(:,2) > h*0.12 & points(:,2) < h*0.88;
    points = points(ids,:);

    [Rc_w, Tc_w] = estimate_cam_pose(sensor, cam_invK);
    [Xc, Yc, Zc, ~] = get_XYZc(Rc_w, Tc_w, points, invK);
end

setPoints(pointTracker, points);
% back up current frame
framei = framei + 1;
oldXc = Xc;
oldYc = Yc;
oldZc = Zc;
oldSensor = sensor;
olddt = dt;

% converting the velocities from the camera frame to the world frame
vel = Rc_w'*v_omega_c(1:3);
omg = Rc_w'*v_omega_c(4:6);

a=1/exp(norm(vel-oldv));
vel = vel*a + oldv*(1-a);
a=1/exp(norm(omg-oldomg)/100);
omg = omg*a + oldomg*(1-a);
oldv = vel;
oldomg = omg;
end
