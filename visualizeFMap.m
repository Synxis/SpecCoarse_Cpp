clc; clear all; close all;
addpath(genpath('~/Dropbox/gpUtils/'))
fileName = 'output_fmap.txt';
fMap = table2array(readtable(fileName));
figure()
plotFMap(fMap)