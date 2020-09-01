function spikeSortingPipeline(inp)
% A function for re-ordering channels and running kilosort.
% spikeSortingPipeline(inp)
% Input: A structure with the following fields:
%        io.dataFiles - raw data file names. If multiple, have to be in a
%          cell array.
%        io.procFolders - kilosort processing folders for storing
%          temp_wh.dat files. If you are stitching input data files
%          together, only a single processing folder is needed. Otherwise,
%          the number of processing folders should correspond to the number
%          of input files. If you do not supply processing folder paths,
%          the output folders will be used as processing folders (the last
%          one in case you are stitching). The intention should be to place
%          processing folders on a local SSD drive, to increase the
%          efficiency of the kilosort algorithm.
%        io.outputFolders - output folders for saving data files with
%          channels re-ordered and/or deleted and possibly medians
%          subtracted depending on the tasks. Each input file has to be
%          matched by an output folder. If you are merging files, there
%          also has to be an extra ouput folder specified for stitched
%          data. Other output files produced by automated Kilosort spike
%          sorting and electrode drift plots will also be saved in
%          corresponding output folders. Folders where Kilosort will run
%          should be free from results from any previous Kilosort runs,
%          which might create problems.
%        io.deleteDataFiles - logical that if true, the files in the
%          dataFiles input cell array will be deleted; false by default.
%        conf.probe - a string specifying the probe type. Currently
%          supported probes are Neuropixels, A32-A1x32-Edge-5mm-20-177,
%          A32-A1x32-5mm-25-177, A32-Buzsaki32-5mm-BUZ-200-160,
%          A32-A1x32-Poly3-5mm-25s-177, A32-A1x32-Poly3-10mm-50-177,
%          A64-Buzsaki64-5mm-BUZ-200-160, A64-A4x4-tet-5mm-150-200-121,
%          CM16LP-A2x2-tet-3mm-150-150-121, CM16LP-A4x4-3mm-100-125-177,
%          CM16LP-A1x16-Poly2-5mm-50s-177, CM16-A1x16-5mm-25-177,
%          CM32-A32-Poly2-5mm-50s-177, CM32-A32-Poly3-5mm-25s-177,
%          CM32-A1x32-6mm-100-177, CM32-A1x32-Edge-5mm-100-177,
%          H32-A1x32-Edge-5mm-20-177, and H32-Buzsaki32-5mm-BUZ-200-160.
%        conf.probeFlip - a logical that is true if the probe/adaptor was
%          connected to the headstage upside-down during the recording
%          session (the labels on the headstage and probe connectors facing
%          opposite sides); default is false.
%        conf.headstage - a string specifying the type of headstage used in
%          combination with the probe. Supported headstages are
%          RHD2132_16ch, RHD2132_32ch, RHD2164_top, RHD2164_bottom,
%          RHD2164, and Neuropixels.
%        conf.nChans - a cell array with the first element being an EEG
%          data channel configuration vector indicating which channels from
%          the original file are contained within the current data file. If
%          full original file is used, then the vector simply corresponds
%          to the original channels (1:end). The second element in the
%          array corresponds to the number of extra input channels that are
%          not electrode recordings. If the cell array is left empty, the
%          default number of EEG recording sites will be assumed based on
%          the probe configuration (conf.probe). If only a single element
%          is supplied, it will be assumed to correspond to the EEG
%          channels only.
%        conf.samplingFrequency (default: 30000).
%        conf.tempFact is the multiplication factor used to determine the
%          number of spike sorting templates to be used by kilosort. The
%          number of templates will be the multiple of conf.tempFact and
%          the number of recording channels. Default is 6.
%        tasks.reorder - a logical that is true for re-ordering channels;
%          true by default.
%        tasks.subtractMedian - a logical that is true for subtracting the
%          median recording trace; false by default
%        tasks.stitchFiles - a logical that if true, input data files are
%          stitched together; false by default. Note that If you are
%          merging files, there has to be an extra ouput folder specified
%          for stitched data.
%        tasks.chanMap - a logical that if true, creates a channel map file
%          (forPRB...); false by default.
%        tasks.runKS - a scalar indicating whether to run automated
%          kilosort and which version. Available options are:
%            1 - kilosort 1;
%            2 - kilosort 2;
%            0 - don't run kilosort (default).
%        tasks.driftPlot - a logical that if true, displays and saves
%          electrode drift plots (relevant only if tasks.runKS was true);
%          false by default.
%        tasks.deleteChans - a vector with channels to be deleted
%          (specified according to the original order).



%% Test user input
if ~isfield(inp, 'io') || isempty(inp.io)
  errMsg = 'data file input and output folders not specified. An example of setting up function input is at the bottom of the file';
  error(['spikeSortingPipeline: ' errMsg])
else
  io = inp.io;
end

if ~isfield(inp, 'conf') || isempty(inp.conf)
  errMsg = 'probe and headstage names not specified. An example of setting up function input is at the bottom of the file';
  error(['spikeSortingPipeline: ' errMsg])
else
  conf = inp.conf;
end

if ~isfield(inp, 'tasks') || isempty(inp.tasks)
  errMsg = 'no tasks specified. An example of setting up function input is at the bottom of the file';
  error(['spikeSortingPipeline: ' errMsg])
else
  tasks = inp.tasks;
  if ~isfield(tasks, 'reorder') && ~isfield(tasks, 'subtractMedian') && ~isfield(tasks, 'stitchFiles') &&...
     ~isfield(tasks, 'chanMap') && ~isfield(tasks, 'runKS') && ~isfield(tasks, 'driftPlot') && ~isfield(tasks, 'deleteChans')
    error('spikeSortingPipeline: no tasks specified')
  end
end


if ~isfield(io, 'dataFiles') || isempty(io.dataFiles)
  error('spikeSortingPipeline: io.dataFiles input not provided')
else
  dataFiles = io.dataFiles;
  if ~iscell(dataFiles)
    dataFiles = {dataFiles};
  end
  for iFile = 1:numel(io.dataFiles)
    if ~(strcmpi(io.dataFiles{iFile}(end-2:end),'dat') || strcmpi(io.dataFiles{iFile}(end-2:end),'bin'))
      errMsg = 'incorrect data file format. Other than dat and bin file formats are not accepted';
      error(['spikeSortingPipeline: ' errMsg])
    end
  end
end

if ~isfield(io, 'outputFolders') || isempty(io.outputFolders)
  error('spikeSortingPipeline: io.outputFolders input not provided')
else
  outputFolders = io.outputFolders;
  if ~iscell(outputFolders)
    outputFolders = {outputFolders};
  end
end
for i1 = 1:numel(outputFolders)-1
  for i2 = 1:i1-1
    assert(~strcmpi(outputFolders{i1}, outputFolders{i2}), ...
      'spikeSortingPipeline: using the same output folder several times creates risk of files being overwritten')
  end
end
clear i1 i2

if ~isfield(io, 'procFolders') || isempty(io.procFolders)
  procFolders = outputFolders;
else
  procFolders = io.procFolders;
  if ~iscell(procFolders)
    procFolders = {procFolders};
  end
end

if ~isfield(io, 'deleteDataFiles') || isempty(io.deleteDataFiles)
  deleteDataFiles = false;
else
  deleteDataFiles = io.deleteDataFiles;
end


if ~isfield(conf, 'probe') || isempty(conf.probe)
  error('spikeSortingPipeline: conf.probe input not provided')
else
  probe = conf.probe;
  if ~strcmpi(probe,'A32-A1x32-Edge-5mm-20-177') && ~strcmpi(probe,'A32-A1x32-5mm-25-177') &&...
     ~strcmpi(probe,'A32-Buzsaki32-5mm-BUZ-200-160') && ~strcmpi(probe,'A32-A1x32-Poly3-5mm-25s-177') &&...
     ~strcmpi(probe,'A32-A1x32-Poly3-10mm-50-177') && ~strcmpi(probe,'A64-Buzsaki64-5mm-BUZ-200-160') &&...
     ~strcmpi(probe,'A64-A4x4-tet-5mm-150-200-121') && ~strcmpi(probe,'CM16LP-A2x2-tet-3mm-150-150-121') &&...
     ~strcmpi(probe,'CM16LP-A4x4-3mm-100-125-177') && ~strcmpi(probe,'CM16LP-A1x16-Poly2-5mm-50s-177') &&...
     ~strcmpi(probe,'CM16-A1x16-5mm-25-177') && ~strcmpi(probe,'CM32-A32-Poly2-5mm-50s-177') &&...
     ~strcmpi(probe,'CM32-A32-Poly3-5mm-25s-177') && ~strcmpi(probe,'CM32-A1x32-6mm-100-177') &&...
     ~strcmpi(probe,'CM32-A1x32-Edge-5mm-100-177') && ~strcmpi(probe,'H32-A1x32-Edge-5mm-20-177') &&...
     ~strcmpi(probe,'H32-Buzsaki32-5mm-BUZ-200-160') && ~strcmpi(probe,'Neuropixels')
    errMsg = ['probe ' probe ' is not supported. Currently supported probes are A32-A1x32-Edge-5mm-20-177, '...
      'A32-A1x32-5mm-25-177, A32-Buzsaki32-5mm-BUZ-200-160, A32-A1x32-Poly3-5mm-25s-177, A32-A1x32-Poly3-10mm-50-177, '...
      'A64-Buzsaki64-5mm-BUZ-200-160, A64-A4x4-tet-5mm-150-200-121, CM16LP-A2x2-tet-3mm-150-150-121, '...
      'CM16LP-A4x4-3mm-100-125-177, CM16LP-A1x16-Poly2-5mm-50s-177, CM16-A1x16-5mm-25-177, CM32-A1x32-6mm-100-177, '...
      'CM32-A32-Poly2-5mm-50s-177, CM32-A32-Poly3-5mm-25s-177, CM32-A1x32-Edge-5mm-100-177, H32-A1x32-Edge-5mm-20-177, '...
      'H32-Buzsaki32-5mm-BUZ-200-160, and Neuropixels'];
    error(['spikeSortingPipeline: ' errMsg])
  end
end

if ~isfield(conf, 'probeFlip') || isempty(conf.probeFlip)
  probeFlip = false;
else
  probeFlip = conf.probeFlip;
end

if ~isfield(conf, 'headstage') || isempty(conf.headstage)
  if strcmpi(probe(1:4), 'CM16')
    headstage = 'RHD2132_16ch';
  elseif strcmpi(probe(1:3), 'A32')
    headstage = 'RHD2164_top';
  elseif strcmpi(probe(1:3), 'A64')
    headstage = 'RHD2164';
  elseif strcmpi(probe, 'Neuropixels')
    headstage = 'Neuropixels';
  else
    error('spikeSortingPipeline: conf.headstage input not provided')
  end
else
  headstage = conf.headstage;
  if ~strcmpi(headstage,'RHD2132_16ch') && ~strcmpi(headstage,'RHD2132_32ch') && ~strcmpi(headstage,'RHD2164_top') &&...
     ~strcmpi(headstage,'RHD2164_bottom') && ~strcmpi(headstage,'RHD2164') && ~strcmpi(headstage,'Neuropixels')
    errMsg = ['headstage ' headstage ' is not supported. Currently supported headstages are RHD2132_16ch, RHD2132_32ch, '...
      'RHD2164_top, RHD2164_bottom, RHD2164, and Neuropixels'];
    error(['spikeSortingPipeline: ' errMsg])
  end
  if (strcmp(probe,'Neuropixels') && ~strcmp(headstage,'Neuropixels')) ||...
      (strcmpi(probe(1:3),'A64') && ~strcmp(headstage,'RHD2164')) ||...
      (strcmpi(probe(1:3),'A32') && (~strcmp(headstage,'RHD2164_top') && ~strcmp(headstage,'RHD2164_bottom') && ~strcmp(headstage,'RHD2164'))) ||...
      ((strcmpi(probe(1:3),'H32') || strcmpi(probe(1:4),'CM32')) &&...
      (~strcmp(headstage,'RHD2132_32ch') && ~strcmp(headstage,'RHD2164_top') && ~strcmp(headstage,'RHD2164_bottom'))) ||...
      (strcmpi(probe(1:4), 'CM16') && ~strcmp(headstage,'RHD2132_16ch'))
    error('spikeSortingPipeline: your probe and headstage are incompatible')
  end
end

if ~isfield(conf, 'nChans') || isempty(conf.nChans)
  chansIgnore = false;
  nChans = [];
else
  chansIgnore = true;
  nChans = conf.nChans;
end

if ~isfield(conf, 'samplingFrequency') || isempty(conf.samplingFrequency)
  samplingFrequency = 30000;
else
  samplingFrequency = conf.samplingFrequency;
end

if ~isfield(conf, 'tempFact') || isempty(conf.tempFact)
  tempFact = 6;
else
  tempFact = conf.tempFact;
end


if ~isfield(tasks, 'reorder') || isempty(tasks.reorder)
  reorder = false;
else
  reorder = tasks.reorder;
end

if ~isfield(tasks, 'subtractMedian') || isempty(tasks.subtractMedian)
  subtractMedian = false;
else
  subtractMedian = tasks.subtractMedian;
end

if  ~isfield(tasks, 'stitchFiles') || isempty(tasks.stitchFiles)
  stitchFiles = false;
else
  stitchFiles = tasks.stitchFiles;
end
if stitchFiles
  procFolders = {procFolders{end}}; %#ok<CCAT1>
end

if  ~isfield(tasks, 'chanMap') || isempty(tasks.chanMap)
  chanMapFile = false;
else
  chanMapFile = tasks.chanMap;
end

if ~isfield(tasks, 'runKS') || isempty(tasks.runKS)
  runKS = 0;
else
  runKS = tasks.runKS;
end

if ~isfield(tasks, 'driftPlot') || isempty(tasks.driftPlot)
  drift = false;
else
  drift = tasks.driftPlot;
end
if isfield(tasks, 'driftPlot') && tasks.driftPlot && isempty(which('ksDriftmap'))
  addpath(genpath('R:\CSN\Shared\Dynamics\Code\github_cortex-lab_spikes'))
  rmpath(genpath('R:\CSN\Shared\Dynamics\Code\github_cortex-lab_spikes\.git'))
elseif ~isfield(tasks, 'driftPlot')
  tasks.driftPlot = false;
end
if tasks.driftPlot
  assert(isempty(strfind(which('findpeaks'), 'chronux')), 'remove chronux from path, otherwise matlab''s own findpeaks function inaccessible') %#ok<STREMP>
end

if ~isfield(tasks, 'deleteChans') || isempty(tasks.deleteChans)
  deleteChans = [];
else
  deleteChans = tasks.deleteChans;
end

if ~reorder && ~subtractMedian && ~stitchFiles && ~runKS && ~drift && isempty(deleteChans)
  error('spikeSortingPipeline: no tasks specified')
end


%% Re-order channels and run kilosort
try
  reset(gpuDevice);
catch
  %do nothing, keep going
end
dataFilesToDelete = dataFiles;
if stitchFiles
  channelMedian = [];
  medianTrace = [];
  for iFile = 1:numel(dataFiles)
    if ~exist(outputFolders{iFile},'dir')
      mkdir(outputFolders{iFile});
    end
    if chansIgnore
%       [channelMedianFile, medianTraceFile, ~, ops.chanMap, dataFilesFull{iFile}, dataFiles{iFile}, ops.NchanTOT, swapOrder,...
%         probe2headstageConf] = swapCAR(dataFiles{iFile}, probe, probeFlip, headstage, reorder, subtractMedian, chanMapFile,...
%         outputFolders{iFile}, deleteChans, nChans); %#ok<AGROW>
      [channelMedianFile, medianTraceFile, ops.NchanTOT, ops.chanMap, dataFiles{iFile}, ~, ~, swapOrder, probe2headstageConf] = swapCAR(...
        dataFiles{iFile}, probe, probeFlip, headstage, reorder, subtractMedian, chanMapFile, outputFolders{iFile}, deleteChans, nChans); %#ok<*ASGLU,*AGROW>
%       connected = zeros(1,nChans);
%       connected(1:ops.NchanTOT) = ones(1,ops.NchanTOT); %#ok<*NASGU>
    else
      [channelMedianFile, medianTraceFile, ~, ops.chanMap, dataFiles{iFile}, ~, ops.NchanTOT, swapOrder, probe2headstageConf] = swapCAR(...
        dataFiles{iFile}, probe, probeFlip, headstage, reorder, subtractMedian, chanMapFile, outputFolders{iFile}, deleteChans);
    end
    channelMedian = [channelMedian channelMedianFile];
    medianTrace = [medianTrace medianTraceFile];
  end
  if ~exist(outputFolders{iFile+1},'dir')
    mkdir(outputFolders{iFile+1});
  end
  stitchedDataFile = [outputFolders{iFile+1} filesep 'stitchedData.dat'];
  save([stitchedDataFile(1:end-4) '_medianTrace.mat'], 'channelMedian', 'medianTrace', 'swapOrder', 'probe2headstageConf', '-v7.3');
  load(ops.chanMap);
  if chanMapFile
    probeConfFile = [outputFolders{iFile+1} filesep 'forPRB_' probe2headstageConf.probeConf.probe '.mat'];
    save(probeConfFile, 'chanMap', 'chanMap0ind', 'connected', 'shankInd', 'xcoords', 'ycoords', '-v7.3'); %#ok<*USENS>
  else
    probeConfFile = dir([outputFolders{iFile+1} filesep 'forPRB*.mat']);
    probeConfFile = probeConfFile.name;
  end
  ops.chanMap = probeConfFile;
  fix_dat_stitch(dataFiles, ops.NchanTOT, samplingFrequency, stitchedDataFile);
  if runKS
    if runKS == 1
      addpath(genpath('R:\Neuropix\Shared\Code\github_cortex-lab_KiloSort'))
      rmpath(genpath('R:\Neuropix\Shared\Code\github_cortex-lab_KiloSort\.git'))
      rmpath(genpath('R:\Neuropix\Shared\Code\github_cortex-lab_KiloSort\CUDA'))
      if verLessThan('matlab', '9.4') % add CUDA8 mex file directory
        addpath('R:\Neuropix\Shared\Code\github_cortex-lab_KiloSort\CUDA\CUDA8')
      elseif  verLessThan('matlab', '9.6') % add CUDA9 mex file directory
        addpath('R:\Neuropix\Shared\Code\github_cortex-lab_KiloSort\CUDA\CUDA9')
      elseif  verLessThan('matlab', '9.7') % add CUDA10.0 mex file directory
        addpath('R:\Neuropix\Shared\Code\github_cortex-lab_KiloSort\CUDA\CUDA10')
      elseif  verLessThan('matlab', '9.8') % add CUDA10.1 mex file directory
        addpath('R:\Neuropix\Shared\Code\github_cortex-lab_KiloSort\CUDA\CUDA101')
      end
    elseif runKS == 2
      addpath(genpath('R:\Neuropix\Shared\Code\github_cortex-lab_KiloSort2'))
      rmpath(genpath('R:\Neuropix\Shared\Code\github_cortex-lab_KiloSort2\.git'))
      rmpath(genpath('R:\Neuropix\Shared\Code\github_cortex-lab_KiloSort2\CUDA'))
      if verLessThan('matlab', '9.4') % add CUDA8 mex file directory
        addpath('R:\Neuropix\Shared\Code\github_cortex-lab_KiloSort2\CUDA\CUDA8')
      elseif  verLessThan('matlab', '9.6') % add CUDA9 mex file directory
        addpath('R:\Neuropix\Shared\Code\github_cortex-lab_KiloSort2\CUDA\CUDA9')
      elseif  verLessThan('matlab', '9.7') % add CUDA10.0 mex file directory
        addpath('R:\Neuropix\Shared\Code\github_cortex-lab_KiloSort2\CUDA\CUDA10')
      elseif  verLessThan('matlab', '9.8') % add CUDA10.1 mex file directory
        addpath('R:\Neuropix\Shared\Code\github_cortex-lab_KiloSort2\CUDA\CUDA101')
      end
    end
    
    load(ops.chanMap);
    if runKS == 1
      [ops.Nfilt, ops.root, ops.fproc, ops.fbinary] = initKS(ops.NchanTOT, outputFolders{iFile+1},...
        procFolders{1}, stitchedDataFile, true, tempFact);
      try
        ks_master_file;
      catch me
        if exist(ops.fproc, 'file')
          delete(ops.fproc);
        end
        disp(getReport(me))
        throw(me);
      end
    elseif runKS == 2
      [~, ops.root, ops.fproc, ops.fbinary] = initKS(ops.NchanTOT, outputFolders{iFile+1},...
        procFolders{1}, stitchedDataFile, true, tempFact);
      ops.nfilt_factor = tempFact;
      try
        ks_master_file2;
      catch me
        if exist(ops.fproc, 'file')
          delete(ops.fproc);
        end
        disp(getReport(me))
        throw(me);
      end
    end
    fclose('all');
  end
else
  for iFile = 1:numel(dataFiles) %#ok<*UNRCH>
    if ~exist(outputFolders{iFile},'dir')
      mkdir(outputFolders{iFile});
    end
    if chansIgnore
%       [~, ~, ~, ops.chanMap, dataFilesFull{iFile}, dataFiles{iFile}, ops.NchanTOT] = swapCAR(dataFiles{iFile}, probe,...
%         probeFlip, headstage, reorder, subtractMedian, chanMapFile, outputFolders{iFile}, deleteChans, nChans); %#ok<AGROW>
      [~, ~, ops.NchanTOT, ops.chanMap, dataFiles{iFile}] = swapCAR(dataFiles{iFile}, probe, probeFlip, headstage, reorder,...
        subtractMedian, chanMapFile, outputFolders{iFile}, deleteChans, nChans);
%       connected = zeros(1,nChans);
%       connected(1:ops.NchanTOT) = ones(1,ops.NchanTOT);
    else
      [~, ~, ~, ops.chanMap, dataFiles{iFile}, ~, ops.NchanTOT] = swapCAR(dataFiles{iFile}, probe, probeFlip, headstage, reorder,...
        subtractMedian, chanMapFile, outputFolders{iFile}, deleteChans);
    end
    if runKS
      if runKS == 1
        addpath(genpath('R:\Neuropix\Shared\Code\github_cortex-lab_KiloSort'))
        rmpath(genpath('R:\Neuropix\Shared\Code\github_cortex-lab_KiloSort\.git'))
        rmpath(genpath('R:\Neuropix\Shared\Code\github_cortex-lab_KiloSort\CUDA'))
        if verLessThan('matlab', '9.4') % add CUDA8 mex file directory
          addpath('R:\Neuropix\Shared\Code\github_cortex-lab_KiloSort\CUDA\CUDA8')
        elseif  verLessThan('matlab', '9.6') % add CUDA9 mex file directory
          addpath('R:\Neuropix\Shared\Code\github_cortex-lab_KiloSort\CUDA\CUDA9')
        elseif  verLessThan('matlab', '9.7') % add CUDA10.0 mex file directory
          addpath('R:\Neuropix\Shared\Code\github_cortex-lab_KiloSort\CUDA\CUDA10')
        elseif  verLessThan('matlab', '9.8') % add CUDA10.1 mex file directory
          addpath('R:\Neuropix\Shared\Code\github_cortex-lab_KiloSort\CUDA\CUDA101')
        end
      elseif runKS == 2
        addpath(genpath('R:\Neuropix\Shared\Code\github_cortex-lab_KiloSort2'))
        rmpath(genpath('R:\Neuropix\Shared\Code\github_cortex-lab_KiloSort2\.git'))
        rmpath(genpath('R:\Neuropix\Shared\Code\github_cortex-lab_KiloSort2\CUDA'))
        if verLessThan('matlab', '9.4') % add CUDA8 mex file directory
          addpath('R:\Neuropix\Shared\Code\github_cortex-lab_KiloSort2\CUDA\CUDA8')
        elseif  verLessThan('matlab', '9.6') % add CUDA9 mex file directory
          addpath('R:\Neuropix\Shared\Code\github_cortex-lab_KiloSort2\CUDA\CUDA9')
        elseif  verLessThan('matlab', '9.7') % add CUDA10.0 mex file directory
          addpath('R:\Neuropix\Shared\Code\github_cortex-lab_KiloSort2\CUDA\CUDA10')
        elseif  verLessThan('matlab', '9.8') % add CUDA10.1 mex file directory
          addpath('R:\Neuropix\Shared\Code\github_cortex-lab_KiloSort2\CUDA\CUDA101')
        end
      end
      
      if ~chanMapFile
        ops.chanMap = dir([fileparts(dataFiles{iFile}) filesep 'forPRB*.mat']);
        ops.chanMap = ops.chanMap.name;
      end
      load(ops.chanMap);
      if runKS == 1
        [ops.Nfilt, ops.root, ops.fproc, ops.fbinary] = initKS(ops.NchanTOT, outputFolders{iFile},...
          procFolders{iFile}, dataFiles{iFile}, false, tempFact);
        try
          ks_master_file;
        catch me
          if exist(ops.fproc, 'file')
            delete(ops.fproc);
          end
          disp(getReport(me))
          throw(me);
        end
      elseif runKS == 2
        [~, ops.root, ops.fproc, ops.fbinary] = initKS(ops.NchanTOT, outputFolders{iFile},...
          procFolders{iFile}, dataFiles{iFile}, false, tempFact);
        ops.nfilt_factor = tempFact;
        try
          ks_master_file2;
        catch me
          if exist(ops.fproc, 'file')
            delete(ops.fproc);
          end
          disp(getReport(me))
          throw(me);
        end
      end
      fclose('all');
    end
  end
end

fclose('all');
recycleState = recycle('on');
if deleteDataFiles
  for iFile = 1:numel(dataFilesToDelete)
    delete(dataFilesToDelete{iFile});
  end
end
% if chansIgnore
%   for iFile = 1:numel(dataFiles)
%     delete(dataFiles{iFile});
%   end
%   if stitchFiles
%     delete(stitchedDataFile);
%     stitchedDataFile = [outputFolders{iFile+1} filesep 'stitchedDataFull.dat'];
%     fix_dat_stitch(dataFilesFull, nChans, samplingFrequency, stitchedDataFile);
%     for iFile = 1:numel(dataFilesFull)
%       delete(dataFilesFull{iFile});
%     end
%   end
%   % elseif stitchFiles
%   %   for iFile = 1:numel(dataFiles)
%   %     delete(dataFiles{iFile});
%   %   end
% end
recycle(recycleState);
fclose('all');


%% Inspect electrode drift
if drift
  close all
  if stitchFiles
    [spikeTimes, spikeAmps, spikeDepths] = ksDriftmap(outputFolders{end});
    plotDriftmap(spikeTimes, spikeAmps, spikeDepths);
    f1 = gcf;
    hgsave(f1, [outputFolders{end} filesep 'Drift_plot_all_spikes']);
    close(f1);
    plotDriftmap(spikeTimes, spikeAmps, spikeDepths, 'show');
    f2 = gcf;
    hgsave(f2, [outputFolders{end} filesep 'Drift_plot_large_spikes']);
    close(f2);
  else
    for iFile = 1:numel(dataFiles)
      [spikeTimes, spikeAmps, spikeDepths] = ksDriftmap(outputFolders{iFile});
      plotDriftmap(spikeTimes, spikeAmps, spikeDepths);
      f1 = gcf;
      hgsave(f1, [outputFolders{iFile} filesep 'Drift_plot_all_spikes']);
      close(f1);
      plotDriftmap(spikeTimes, spikeAmps, spikeDepths, 'show');
      f2 = gcf;
      hgsave(f2, [outputFolders{iFile} filesep 'Drift_plot_large_spikes']);
      close(f2);
    end
  end
end
return


%% Functions
function [Nfilt, root, fproc, fbinary] = initKS(Nchan, root, proc, dataFile, stitchFiles, tempFact)

if stitchFiles
  Nfilt = ceil(tempFact*(Nchan-1)/32)*32; % number of filters to use (2-4 times more than Nchan, should be divisible by 32)
else
  Nfilt = ceil(tempFact*(Nchan-1)/32)*32; %Nfilt = ceil(3*Nchan/32)*32; % number of filters to use (2-4 times more than Nchan, should be divisible by 32)
end
fproc = fullfile(proc, 'temp_wh.dat');
fbinary = dataFile;
return

%% Example:



inp.io.dataFiles = {'R:\Neuropix\Shared\Data\M191018_MD\original_data\TCB-2_g0_t0.imec0.ap.bin'};
inp.io.procFolders = {'D:\'};
inp.io.outputFolders = {'R:\Neuropix\Shared\Data\M191018_MD'};
inp.io.deleteDataFiles = false;
inp.conf.probe = 'Neuropixels';
inp.conf.probeFlip = false;
inp.conf.headstage = 'Neuropixels';
inp.conf.nChans = {1:384; 1};
inp.conf.samplingFrequency = 30000;
inp.conf.tempFact = 3;
inp.tasks.reorder = false;
inp.tasks.subtractMedian = false;
inp.tasks.stitchFiles = false;
inp.tasks.chanMap = true;
inp.tasks.runKS = 1;
inp.tasks.driftPlot = true;
inp.tasks.deleteChans = [];

spikeSortingPipeline(inp);
