% default options are in parenthesis after the comment

try
  gpuArray(1);
  ops.GPU               = 1; % whether to run this code on an Nvidia GPU (much faster, mexGPUall first)
catch
  ops.GPU               = 0;
end
ops.parfor              = 0; % whether to use parfor to accelerate some parts of the algorithm
ops.verbose             = 1; % whether to print command line progress
ops.showfigures         = 1; % whether to plot figures during optimization

if strcmpi(probe, 'Neuropixels')
  ops.datatype            = 'bin';  % binary ('dat', 'bin') or 'openEphys'
else
  ops.datatype            = 'dat';
end

ops.fs                  = 30000;        % sampling rate
% ops.NchanTOT            = 32;           % total number of channels
% ops.Nchan               = sum(connected);           % number of active channels 
% ops.Nfilt               = 768;           % number of filters to use (2-4 times more than Nchan, should be a multiple of 32)     
ops.nNeighPC            = 4; % visualization only (Phy): number of channnels to mask the PCs, leave empty to skip (12)
ops.nNeigh              = 16; % visualization only (Phy): number of neighboring templates to retain projections of (16)

% options for channel whitening
ops.whitening           = 'full'; % type of whitening (default 'full', for 'noSpikes' set options for spike detection below)
ops.nSkipCov            = 1; % compute whitening matrix from every N-th batch (1)
ops.whiteningRange      = 32; % how many channels to whiten together (Inf for whole probe whitening, should be fine if Nchan<=32)

% define the channel map as a filename (string) or simply an array
% ops.chanMap             = 'C:\DATA\Spikes\Piroska\chanMap.mat'; % make this file using createChannelMapFile.m
ops.criterionNoiseChannels = 0.2; % fraction of "noise" templates allowed to span all channel groups (see createChannelMapFile for more info).
% ops.chanMap = 1:ops.Nchan; % treated as linear probe if a chanMap file

% other options for controlling the model and optimization
ops.Nrank               = 3;    % matrix rank of spike template model (3)
ops.nfullpasses         = 6;    % number of complete passes through data during optimization (6)
ops.maxFR               = 20000;  % maximum number of spikes to extract per batch (20000)
ops.fshigh              = 300;   % frequency for high pass filtering
ops.ntbuff              = 64;    % samples of symmetrical buffer for whitening and spike detection
ops.scaleproc           = 200;   % int16 scaling of whitened data
ops.NT                  = 32*1024+ ops.ntbuff;% this is the batch size (try decreasing if out of memory) 
% for GPU should be multiple of 32 + ntbuff

% the following options can improve/deteriorate results. 
% when multiple values are provided for an option, the first two are beginning and ending anneal values, 
% the third is the value used in the final pass. 
ops.Th               = [4 10 10];    % threshold for detecting spikes on template-filtered data ([6 12 12])
ops.lam              = [5 20 20];   % large means amplitudes are forced around the mean ([10 30 30])
ops.nannealpasses    = 4;            % should be less than nfullpasses (4)
ops.momentum         = 1./[20 400];  % start with high momentum and anneal (1./[20 1000])
ops.shuffle_clusters = 1;            % allow merges and splits during optimization (1)
ops.mergeT           = .1;           % upper threshold for merging (.1)
ops.splitT           = .1;           % lower threshold for splitting (.1)

% options for initializing spikes from data
ops.initialize      = 'no'; %'fromData' or 'no'
ops.spkTh           = -6;      % spike threshold in standard deviations (4)
ops.loc_range       = [3  1];  % ranges to detect peaks; plus/minus in time and channel ([3 1])
ops.long_range      = [30  6]; % ranges to detect isolated peaks ([30 6])
ops.maskMaxChannels = 5;       % how many channels to mask up/down ([5])
ops.crit            = .65;     % upper criterion for discarding spike repeates (0.65)
ops.nFiltMax        = 10000;   % maximum "unique" spikes to consider (10000)

% load predefined principal components (visualization only (Phy): used for features)
dd                  = load('PCspikes2.mat'); % you might want to recompute this from your own data
ops.wPCA            = dd.Wi(:,1:7);   % PCs 

% options for posthoc merges (under construction)
ops.fracse  = 0.1; % binning step along discriminant axis for posthoc merges (in units of sd)
ops.epu     = Inf;

ops.ForceMaxRAMforDat   = 20e9; % maximum RAM the algorithm will try to use; on Windows it will autodetect.

%% 
tic; % start timer

if strcmp(ops.datatype , 'openEphys')
   ops = convertOpenEphysToRawBInary(ops);  % convert data, only for OpenEphys
end
%
[rez, DATA, uproj] = preprocessData(ops);

if strcmp(ops.initialize, 'fromData')
    % do scaled kmeans to initialize the algorithm (not sure if functional yet for CPU)
    optimizePeaks(uproj);
end
%
rez = fitTemplates(rez, DATA, uproj);

%
% extracts final spike times (overlapping extraction)
gpuDevice(1);
rez                = fullMPMU(rez, DATA);

% posthoc merge templates (under construction)
%     rez = merge_posthoc2(rez);

% save matlab results file
save(fullfile(ops.root,  'rez.mat'), 'rez', '-v7.3');

% save python results file for Phy
rezToPhy(rez, ops.root);

% remove temporary file
delete(ops.fproc);
%%
