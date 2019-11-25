% Demo script for IDvortex

% Prepare:
clc; clear; close all; tic;
MatLabSettings

% Load sample data:
load('sampleXYUV.mat')
% x,y coordinate vectors (mm) and Vx and Vy fields (m/s)

% Call the IDvortex scripts
vortex=IDvortex(x,y,U,V);

% Plot the U field with the found vortex
figure; surface('ZData',U,'YData',y,'XData',x,'CData',U,'FaceColor','interp','EdgeColor','none'); title('U')
hold on; plot3(vortex(:,1),vortex(:,2),9e9,'g+','MarkerSize',20)

% Plot the V field with the found vortex
figure; surface('ZData',V,'YData',y,'XData',x,'CData',V,'FaceColor','interp','EdgeColor','none'); title('V')
hold on; plot3(vortex(:,1),vortex(:,2),9e9,'g+','MarkerSize',20)

toc