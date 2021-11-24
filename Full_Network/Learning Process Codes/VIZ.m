%% Vizualiation Code
clear all; close all; clc;

%% Open File and Ploat

DATA = load('UNIFULL.txt');

X = DATA(:,1); Y = DATA(:,2);

plot(X,Y)