clc
clear all
close all
format longG
Title ={'x轴','y轴','z轴','x轴','y轴','z轴'};
Label ={'位置误差（m）','位置误差（m）','位置误差（m）','速度误差（m/s）','速度误差（m/s）','速度误差（m/s）'};

% Add path to the Earth plotting function. 
addpath([pwd, '/PlotEarth/PlotEarth']);

% Add path to the HPOP-propagator. 
addpath([pwd, '/HPOP-propagator']);

% Add path to the ukf. 
addpath([pwd, '/ukf']);

% Add path to the TLE data (these are just some examples with a focus on
% GNSS). 
%addpath([pwd, '/TLE_Files']);

%% LOAD PHYSICAL CONSTANTS INTO THE GLOBAL WORKSPACE

physical_constants

global mu

%% SELECT THE TLE FILES TO PLOT

%TEST
filenames = {'TEST'};


colors = lines(length(filenames));
% colors = hsv(length(filenames));

%% Plot The Earth

% Plot the Earth. 
% If you want a color Earth, use 'neompa', 'BlueMarble'.
% If you want a black and white Earth, use 'neomap', 'BlueMarble_bw'.
% A smaller sample step gives a finer resolution Earth.
h = plotearth('neomap', 'BlueMarble', 'SampleStep', 2);



%% load the real state and compute mearsurement
step=1;

D_real=load('D_real.txt');
D_real=D_real(2:step:end,:);
%D_real(:,4:6)=D_real(:,4:6)*1e3;
 %D_real=load('dst1.mat');
 %D_real=D_real.D_real;
 %D_real=D_real(1:end-1,2:end);
f = 1500;
dx = 1;
res = 512;
name = dir('E:\实验室\高轨定姿\screen\screen\3\*.bmp');
path = 'E:\实验室\高轨定姿\screen\screen\3\';
rect = [];
img_pos = zeros(2,size(name,2));
Z = zeros(3,size(name,2));
figure(2);
hold on
for i=1:size(name,1)
img = imread(strcat(path,name(i).name));
[img_crop,rect] = imcrop(img,rect);
[tic_x,tic_y] = find_center(img_crop);
img_pos(:,i) = [tic_x,tic_y]'+[rect(1),rect(2)]'-1;
Z(:,i) = [(img_pos(1,i)-res/2)*dx;(img_pos(2,i)-res/2)*dx;f];
Z(:,i) = Z(:,i)./norm(Z(:,i));
imshow(img);
hold on
plot(img_pos(1,i),img_pos(2,i),'r.'); 
pause(0.3)
rect(1,1:2) = round(img_pos(:,i)')-[rect(3),rect(4)]/2;
end



S_real=load('sat.txt');
S_real=S_real(2:step:end,:);
%mearsurement
Z = (D_real(:,1:3)-S_real(:,1:3));
%normlize
Z = Z./(sqrt(Z(:,1).^2+Z(:,2).^2+Z(:,3).^2)*[1 1 1]);
%add noise
Z_noise = Z + 0.0001*rand(size(Z,1),3);



%plot
plot3(D_real(:,1) / R_e, D_real(:,2) / R_e, D_real(:,3) / R_e,...
            'color', colors(1,:), 'LineWidth', 1)
        plot3(D_real(1,1) / R_e, D_real(1,2) / R_e, D_real(1,3) / R_e,...
            '.', 'color', colors(1,:), 'MarkerSize', 10)
plot3(S_real(:,1) / R_e, S_real(:,2) / R_e, S_real(:,3) / R_e,...
            'color', colors(1,:), 'LineWidth', 1)
        plot3(S_real(1,1) / R_e, S_real(1,2) / R_e, S_real(1,3) / R_e,...
            '.', 'color', colors(1,:), 'MarkerSize', 10)
        hold on

%% Determinant The Inial Orbit
%start time is 1
S0 = Gauss(Z([1 3 5],:),S_real([1 3 5],1:3),[0,2*step,4*step]);
%initialize the propagtor
init
%% UKF Determinantion 
P = diag([1e5 1e5 1e5 0.01 0.01 0.001]);
Q = diag([0.1 0.1 0.1 0.001 0.001 0.001]);
R = diag([0.0001 0.0001 0.0001]);
f_func = @propagator;
h_func = @observation;
M = S0';
U_MM = zeros(size(M,1),size(Z,1));
U_PP = zeros(size(M,1),size(M,1),size(Z,1));
handle=waitbar(0,'wait...');
U_MM(:,3)   = M;
U_PP(:,:,3) = P;
for k=4:size(Z,1)
   t0=tic;
   [M,P,X_s,w] = ukf_predict3(M,P,f_func,Q,R*eye(1),step,0.0001,2,[],0);
   [M,P] = ukf_update3(M,P,Z(k,:)',h_func,R*eye(1),X_s,w,S_real(k,:),0.0001,2,[],0);
   U_MM(:,k)   = M;
   U_PP(:,:,k) = P;
   t=toc(t0);
   t=(size(Z,1)-k)*t;
   waitbar(k/size(Z,1),handle,sprintf('%.2f%%,Remaining time : %.1fmin%.1fs',k/size(Z,1)*100,floor(t/60),mod(t,60)));
end
err=D_real'-U_MM;
for i=1:6
figure(i+1)
plot((1:(size(err,2)-2))*step,err(i,3:end));
set(0,'defaultfigurecolor','w');
box off;
title(Title(i));
xlabel('时间（s）');
ylabel(Label(i));
end

close(handle);

%[U_SM, U_SP] = urts_smooth1(U_MM,U_PP,f_func,Q,dt);





