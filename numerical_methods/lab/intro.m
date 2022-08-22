clc
clear all
close all

iex = 1; % Select example
%% EXAMPLE 1: Projectile
if iex==1
    t = linspace(0,100,5);   % in seconds
    v0 = 500;                   % in m/s
    theta = 30;                 % in degrees
    
    % constant
    g = 9.8;                    % m/s2 (gravity constant)
    %===============================================
      
    % horizontal and vertical displacements  % KEY
    x = v0 * t * cosd(theta); % 500 x [1 2 3 4.... 100] * cosd(30)
    y = v0 * t * sind(theta) - 0.5 * g * t.^2;
    
    %===============================================
    % Plot
    figure
    plot(x,y,'LineWidth',3)
    xlabel('x'), ylabel('y')
    
    % LAB ASSIGNMENT 1
    % Note that x, y has three variable (v0, t, theta)
    % 1. Make inline functions for x and y
    % 2. Plot for theta = 30, 45, 80 degrees    
end

%% EXAMPLE 2: Read 1D waveform data
if iex == 2
    %ddir = '~/REPOSITORIES/IITR_seismo/classes/ES510/lab/data/';
    ddir = '/Users/vipul/Desktop/classes/numerical_methods/lab/data/';
    fname = 'XZ_GOAT_BHZ.dat';
    fid = fopen(strcat(ddir,fname));
    
    % OPTION 1
    data = textscan(fid,'%f');
    data_arr = data{1};     
    % OPTION 2
    data2 = dlmread(strcat(ddir,fname));
    
    % Plot
    figure
    plot(data_arr)
    figure
    plot(data2)
    
    % Add information from the header
    % Standard seismic waveform readers (obspy, SAC) will read this information
    % direclty
    
    % LAB ASSIGNMENT 2
    % Find x where y is minimum and maximum
end
