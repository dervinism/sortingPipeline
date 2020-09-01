%script for running Kilosort

%%
% mouseName = 'ALK052';
% thisDate = datestr(now, 'yyyy-mm-dd');
% thisDate = '2017-08-23';
% fnBase = [mouseName '_' thisDate '_g0_t0.imec.'];
fnBase = 'continuous_swappedNoCAR';
%% parameters

%ops.chanMap = 'forPRBimecP3opt3.mat';
ops.chanMap = 'forPRB_A2x2_tet_3mm_150_150_121.mat';
ops.NchanTOT = 16;
ops.Nfilt = 4*16; % number of filters to use (2-4 times more than Nchan, should be a multiple of 32)     
% ops.root = ['D:\Data', filesep, mouseName, filesep, thisDate, filesep];
ops.root = 'R:\CSN\Shared\Dynamics\spikeSorting\';
ops.fproc = fullfile('R:\CSN\Shared\Dynamics\spikeSorting', 'temp_wh.dat');

% fn = fullfile([ops.root fnBase '.bin']);
% fnAfterCAR = fullfile([ops.root fnBase '.bin']);
fn = fullfile([ops.root fnBase '.dat']);
fnAfterCAR = fullfile([ops.root fnBase '.dat']);

load(ops.chanMap);

ext = fn(end-2:end);

ops.fbinary = fnAfterCAR;


% % first perform CAR
% tic
% medianTrace = applyCARtoDat(fn, ops.NchanTOT, ['/localdisk/mush/' mouseName filesep thisDate]);
% toc

%% then run KS
ks_master_file;
fclose('all');
return

% %% then copy to server
% tic
% % basketDrive = 'Z:\';
% % zserverDrive = 'X:\';
% % lugaroDrive = 'Y:\';
% basketDrive = '\\basket.cortexlab.net\data\';
% zserverDrive = '\\zserver.cortexlab.net\data\';
% lugaroDrive = '\\lugaro.cortexlab.net\staging\';
% lugaroDrive2 = '\\lugaro.cortexlab.net\toarchive\';
% 
% zserverDest = fullfile(zserverDrive, 'multichanspikes', mouseName, thisDate);
% basketDest = fullfile(basketDrive, 'nick', mouseName, thisDate);
% lugaroDest = fullfile(lugaroDrive, [mouseName '_' thisDate]);
% 
% % npy files go to basket
% mkdir(basketDest)
% fprintf(1, 'moving npy to basket\n');
% movefile(fullfile(ops.root, '*.npy'), basketDest);
% movefile(fullfile(ops.root, 'params.py'), basketDest);
% movefile(fullfile(ops.root, 'rez.mat'), basketDest);
% toc
% %% afterCAR (along with LFP, meta, median) goes to zserver
% mkdir(zserverDest);
% fprintf(1, 'moving CAR to zserver\n');
% movefile(fnAfterCAR, zserverDest);
% toc
% 
% %%
% fprintf(1, 'moving other to zserver\n');
% movefile(fullfile(ops.root, [fnBase 'ap.meta']), zserverDest);
% movefile(fullfile(ops.root, [fnBase 'ap_medianTrace.mat']), zserverDest);
% movefile(fullfile(ops.root, [fnBase 'lf.bin']), zserverDest);
% movefile(fullfile(ops.root, [fnBase 'lf.meta']), zserverDest);
% toc
% %% raw goes to toarchive
% fprintf(1, 'moving raw to lugaro\n');
% tic
% mkdir(lugaroDest);
% movefile(fullfile(ops.root, [fnBase 'ap.bin']), lugaroDest);
% toc
% fprintf(1, 'copying to lugaro''s toarchive\n');
% tic
% movefile(lugaroDest, lugaroDrive2);
% toc
% fprintf(1, 'done!\n');
% toc