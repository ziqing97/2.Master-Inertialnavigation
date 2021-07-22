clc
close all
clear all

%%
vn310gnss = importgnss("E:\Studium\M2-Inertialnavigation\Uebung\3\Code\vn310-gnss.csv", [2, Inf]);
vn310imu = importimu("E:\Studium\M2-Inertialnavigation\Uebung\3\Code\vn310-imu.csv", [2, Inf]);