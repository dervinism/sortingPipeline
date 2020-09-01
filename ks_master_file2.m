ops.GPU                = 1; % whether to run this code on an Nvidia GPU (much faster, mexGPUall first)
ops.fs                 = 30000;        % sampling rate
ops.nSkipCov           = 25; % compute whitening matrix from every N-th batch (1)
ops.whiteningRange     = 32; % how many channels to whiten together (Inf for whole probe whitening, should be fine if Nchan<=32)
ops.fshigh             = 150;   % frequency for high pass filtering
ops.ntbuff             = 64;    % samples of symmetrical buffer for whitening and spike detection
ops.scaleproc          = 200;   % int16 scaling of whitened data
ops.NT                 = 64*1024+ ops.ntbuff;% this is the batch size (try decreasing if out of memory)
ops.Th                 = [10 2];        % threshold for detecting spikes on template-filtered data ([6 12 12])
ops.lam                = 10;           % large means amplitudes are forced around the mean ([10 30 30])
ops.momentum           = [20 400];     % start with high momentum and anneal (1./[20 1000])
ops.spkTh              = -6;      % spike threshold in standard deviations (4)
ops.trange             = [0 Inf]; % time range to sort
ops.ThPre              = 8; % threshold crossings for pre-clustering (in PCA projection space)
ops.minfr_goodchannels = 0.1; % minimum firing rate on a "good" channel (0 to skip)
ops.AUCsplit           = 0.9; % splitting a cluster at the end requires at least this much isolation for each sub-cluster (max = 1)
ops.minFR              = 1/50; % minimum spike rate (Hz), if a cluster falls below this for too long it gets removed
ops.sigmaMask          = 30; % spatial constant in um for computing residual variance of spike
ops.nPCs               = feature('numcores'); % how many PCs to project the spikes into
ops.useRAM             = 0; % not yet available

%% 
tic; % start timer

gpuDevice(1);   %re-initialize GPU

rez = preprocessDataSub(ops);

% time-reordering as a function of drift
rez = clusterSingleBatches(rez);

% main tracking and template matching algorithm
rez = learnAndSolve8b(rez);

% final merges
rez = find_merges(rez, 1);

% final splits by SVD
rez = splitAllClusters(rez, 1);

% final splits by amplitudes
rez = splitAllClusters(rez, 0);

% decide on cutoff
rez = set_cutoff(rez);

fprintf('found %d good units \n', sum(rez.good>0))

% write to Phy
fprintf('Saving results to Phy  \n')
rezToPhy(rez, ops.root);

%% if you want to save the results to a Matlab file... 

% discard features in final rez file (too slow to save)
rez.cProj = [];
rez.cProjPC = [];

% save final results as rez2
fprintf('Saving final results in rez2  \n')
fname = fullfile(ops.root, 'rez2.mat');
save(fname, 'rez', '-v7.3');

% remove temporary file
delete(ops.fproc);
%%
