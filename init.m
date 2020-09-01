fclose all;
close all;
clear all %#ok<*CLALL>
clc

inp.io.dataFiles = {'R:\Neuropix\md406\continuous.dat'};
inp.io.procFolders = {'R:\Neuropix\md406'};
inp.io.outputFolders = {'R:\Neuropix\md406'};
inp.io.deleteDataFiles = false;
inp.conf.probe = 'A64-A4x4-tet-5mm-150-200-121';
inp.conf.probeFlip = false;
inp.conf.headstage = 'RHD2164';
inp.conf.nChans = {1:64};
inp.conf.samplingFrequency = 30000;
inp.conf.tempFact = 3;
inp.tasks.reorder = false;
inp.tasks.subtractMedian = true;
inp.tasks.stitchFiles = false;
inp.tasks.chanMap = true;
inp.tasks.runKS = 1;
inp.tasks.driftPlot = true;
inp.tasks.deleteChans = [];

spikeSortingPipeline(inp);