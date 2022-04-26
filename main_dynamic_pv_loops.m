
clc
clear all
close all

%load data
prime_data  = read_PRiME_file(); %load BIN file recording from PRiME
volume_data = loadSuiteHeartVolumeExcel(); %Load excel report with 3D volume generated from long-axis cine using SuiteHeart

data = processPVLoops(prime_data, volume_data);

