%% 
clear;
clc;
addpath(genpath('./'));

%% check data format
load('data/studentdata1.mat');

%%
vels = zeros(3, numel(data));
omgs = zeros(3, numel(data));

vicon_vels = zeros(3, numel(data));
vicon_omgs = zeros(3, numel(data));
ts = zeros(1,numel(data));

pos = zeros(3, numel(data));
ori = zeros(3, numel(data));
vicon_pos = zeros(3, numel(data));
vicon_ori = zeros(3, numel(data));
%%

USE_EKF1 = 0;

tic
% for i = 1:numel(data)
n_tags = zeros(1, numel(data));
for i = 1:numel(data)
    n_tags(i) = numel(data(i).id);
    if mod(i,200) == 0
        disp(int2str(i))
    end
    if numel(data(i).id)
        vicon_v = get_vicon(data(i), vicon, time);
        vic = {};
        vic.t = data(i).t;
        vic.vel = vicon_v(7:12);
%         ekf1 uses vicon's velocities data as measurements
        if USE_EKF1 == 1
            [X, Z, rpy] = ekf1(data(i), vic);
            [vel, omg] = estimate_vel(data(i));
        else 
%         ekf2 uses vison-based velocities estimations as measurements
%         UKF enclosed, please check the flag USE_UKF in ekf2()
            [X, Z, rpy, omg] = ekf2(data(i));
            vel = X(4:6);
        end
        
        vels(:, i) = vel;
        omgs(:, i) = omg;
        pos(:, i) = X(1:3);
        ori(:, i) = rpy;
        vicon_vels(:, i) = vicon_v(7:9);
        vicon_omgs(:, i) = vicon_v(10:end);
        vicon_pos(:,i) = vicon_v(1:3);
        vicon_ori(:,i) = vicon_v(4:6);
    end
    ts(i) = data(i).t;
end
toc

%% 
plot_vel_omg(vels, omgs, vicon_vels, vicon_omgs, ts);
plotpos(pos, vicon_pos, ts, n_tags);
plotrpy(ori, vicon_ori, ts, n_tags);

