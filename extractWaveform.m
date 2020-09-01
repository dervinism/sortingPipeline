function [waveforms, maxWaveforms, cluIDs, spikeCentreIndex, amplitudes, maxChan] = extractWaveform(inp)
% [waveforms, maxWaveforms, cluIDs, spikeCentreIndex, amplitudes, maxChan] = extractWaveform(inp)
% extractWaveform function extracts average waveforms and other associated
% information about spike waveforms following kilosort and phy spike
% sorting.
% Input: inp structure with the following field (some have defaults and are optional)
%        inp.dataDir
%        inp.dataFile
%        inp.chansIgnore - a number of channels to be ignored (not electrode channels), 0 by default
%        inp.outputFile - ('waveforms.mat' by default)
%        inp.display - true or false for displaying average waveforms.
%        inp.wavelength - the duration of the spike waveform (200 by default)
%        inp.merge - merge waveforms for all files: true or false (default).
%                    Set to true to save the results in inp.dataDir
%        inp.dataType ('int16' by default, specify otherwise)
%
% Output: waveforms - a cell array of average waveforms for each stitched
%                     file.
%         maxWaveforms - a cell array of average waveforms containing
%                        channels with the largest amplitude only.
%         cluIDs - a vector of unit IDs corresponding waveforms.
%         spikeCentreIndex - a spike centre index on the waveform.
%         amplitudes - spike amplitudes.
%         maxChan - channel IDs with the largest spike amplitudes.
% The function also saves an output file containing output variables, as
% well as the list of stitched files corresponding to the output data.


%% User input
  dataDir = inp.dataDir; 
  dataFile = inp.dataFile; 

if ~isfield(inp, 'dataType') || isempty(inp.dataType) 
  dataType = 'int16';
else
  dataType = inp.dataType;
end
if ~isfield(inp, 'chansIgnore') || isempty(inp.chansIgnore)
  chansIgnore = 0;
else
  chansIgnore = inp.chansIgnore;
end
if ~isfield(inp, 'outputFile') || isempty(inp.outputFile)
  outputFile = 'waveforms.mat';
else
  outputFile = inp.outputFile;
end
if ~isfield(inp, 'display') || isempty(inp.display)
  display = true;
else
  display = inp.display;
end
if ~isfield(inp, 'wavelength') || isempty(inp.wavelength)
  wavelength = 200;
else
  wavelength = inp.wavelength;
end
if ~isfield(inp, 'merge') || isempty(inp.merge)
  merge = false;
else
  merge = inp.merge;
end


%% Extract the spike cluster info and template waveforms
sp = loadKSdir(dataDir);
cluFull = sp.clu;
resFull = round(sp.st * sp.sample_rate);
clu = [];
res = [];
for sh = 1:16
  if exist([dataDir filesep dataFile(1:end-4) '.clu.' num2str(sh)], 'file') && exist([dataDir filesep dataFile(1:end-4) '.res.' num2str(sh)], 'file')
    resSh = load([dataDir filesep dataFile(1:end-4) '.res.' num2str(sh)]);
    cluSh = load([dataDir filesep dataFile(1:end-4) '.clu.' num2str(sh)]);
    assert(numel(resSh) == numel(cluSh) - 1); cluSh = cluSh(2:end);
    resSh = resSh(cluSh > 0); % removing noise spikes
    cluSh = cluSh(cluSh > 0);
    res = [res; resSh];    
    clu = [clu; cluSh];
  elseif sh == 1
    [clu, res] = resCluFromKilosort(dataDir, 1, 10000, 1:10000);
    clu = clu(2:end);
  end
end
[res, swapOrder] = sort(res);
clu = clu(swapOrder);
resFullDiff = resFull' - [0 resFull(1:end-1)'];
resDiff = res' - [0 res(1:end-1)'];
resDiff = resDiff(2:end);
whichAmps = strfind(resFullDiff,resDiff)-1:strfind(resFullDiff,resDiff)-1+numel(resDiff);
if isempty(whichAmps)
  error('Waveform amplitudes could not be determined');
end
amps = sp.tempScalingAmps(whichAmps);
templateWaveforms = sp.temps;
assert(numel(clu) == numel(res) && numel(res) == numel(amps));
cids = double(sp.cids); cgs = sp.cgs;
cluIDsFull = double(unique(cluFull));
if any(cluIDsFull == 0)
  m = max(cluIDsFull) + 1;
  %clu(clu == 0) = m;
  cids(cids == 0) = m;
  cluIDsFull(cluIDsFull == 0) = m;
end
if any(cluIDsFull == 1)
  m = max(cluIDsFull) + 1;
  %clu(clu == 1) = m;
  cids(cids == 1) = m;
  cluIDsFull(cluIDsFull == 1) = m;
end
assert(max(abs(sort(torow(cluIDsFull)) - sort(torow(cids)))) == 0, 'should be fully compatible');
assert(~any(cgs >= 3), 'unsorted units remain');
cluIDs = double(unique(clu));
cluIDs = cluIDs(cluIDs > 1);
for iClu = 1:numel(cluIDsFull)
  if ~sum(cluIDsFull(iClu) == cluIDs)
    cids(iClu) = NaN;
    cgs(iClu) = NaN;
  end
end
cids = cids(~isnan(cids));
cgs = cgs(~isnan(cgs));
if isempty(cluIDs) && isempty(cids)
  disp('The file contains no single unit activity. Please check with Phy if this is indeed correct. extractWaveform function is terminating.');
  waveforms = []; maxWaveforms = []; cluIDs = []; spikeCentreIndex = []; amplitudes = []; maxChan = [];
  return
end
assert(max(abs(sort(torow(cluIDs)) - sort(torow(cids)))) == 0, 'should be fully compatible');


%% Load raw data
probeConfFile = dir([dataDir filesep 'forPRB*']);
if isempty(probeConfFile.name)
  error('No probe configuration file found in the data folder')
else
  load([dataDir filesep probeConfFile.name],'connected');
end
if exist([dataDir filesep dataFile(1:end-4) '.mat'], 'file')
  load([dataDir filesep dataFile(1:end-4) '.mat']); %#ok<LOAD>
  nFiles = numel(dataPoints); %#ok<NODEF>
end
%sp.n_channels_dat = numel(sp.xcoords)+inp.chansIgnore;

chunkSize = 1000000;
fileName = fullfile(dataDir,dataFile);
fprintf('extractWaveform: working on %s, which is presumed to have %d channels\n', fileName, sp.n_channels_dat)
filenamestruct = dir(fileName);
dataTypeNBytes = numel(typecast(cast(0, dataType), 'uint8')); % The size of a single data point in bytes
nSampsTotal = filenamestruct.bytes/sp.n_channels_dat/dataTypeNBytes;
nChunksTotal = ceil(nSampsTotal/chunkSize);

if ~exist([dataDir filesep dataFile(1:end-4) '.mat'], 'file')
  dataPoints = nSampsTotal;
  nFiles = 1;
  datFileList{1} = fileName;
end

fid = fopen(fileName, 'r');
chunkInd = 1;
%templateLength = size(templateWaveforms,2);
templateLength = wavelength;
spikeCentreIndex = templateLength/2 + 1;
templateCh = sum(connected);
for iFile = 1:nFiles
  templateWaveformsNew{iFile} = nan(numel(cluIDs), templateLength, templateCh); %#ok<*AGROW,*SAGROW>
  nTemplateWaveformsNew{iFile} = nan(numel(cluIDs), templateLength, templateCh);
  amplitudes{iFile} = nan(size(cluIDs));
  nAmplitudes{iFile} = nan(size(cluIDs));
  if merge
    nFiles = 1;
    break %#ok<*UNRCH>
  end
end
while 1
  fprintf(1, 'chunk %d/%d\n', chunkInd, nChunksTotal);
  dat = fread(fid, [sp.n_channels_dat chunkSize], ['*' dataType]);
  if ~isempty(dat)
%     if chansIgnore
%       dat = dat(1:end-chansIgnore,:);
%     end
    
    
    %% Add the spikes to the waveforms
    dataRange = round([(chunkInd-1)*chunkSize+1 chunkInd*chunkSize]);
    %tRange = dataRange./sp.sample_rate;
    resOI = res((res>=dataRange(1)) == (res<=dataRange(2)));
    cluOI = clu((res>=dataRange(1)) == (res<=dataRange(2)));
    ampOI = amps((res>=dataRange(1)) == (res<=dataRange(2)));
    %sum(cluOI == 316)/numel(cluOI)
    for i = 1:numel(resOI)
      iWaveform = zeros(1,templateLength);
      waveform = zeros(templateCh,templateLength);
      iClu = cluOI(i);
      iCluIDs = find(cluIDs == iClu);
      spikeRangeInit = [resOI(i)-templateLength/2 resOI(i)+templateLength/2-1];
      spikeRange = round(spikeRangeInit);
      for iFile = 1:nFiles
        if iFile == 1
          dataPointStart = 1;
          dataPointEnd = sum(dataPoints);
        else
          dataPointStart = sum(dataPoints(1:iFile-1)) + 1;
          dataPointEnd = sum(dataPoints(1:iFile));
        end
        if resOI(i) <= dataPointEnd
          break
        end
      end
      if merge
        filePointer = 1;
      else
        filePointer = iFile;
      end
      if spikeRange(1) < dataPointStart
        startWaveform = round(dataPointStart - spikeRange(1) + 1);
        spikeRange(1) = round(dataPointStart);
      elseif spikeRange(1) < dataRange(1)
        startWaveform = round(dataRange(1) - spikeRange(1) + 1);
        spikeRange(1) = round(dataRange(1));
      else
        startWaveform = 1;
      end
      if spikeRange(2) > dataPointEnd
        endWaveform = round(templateLength - (spikeRange(2) - dataPointEnd));
        spikeRange(2) = round(dataPointEnd);
      elseif spikeRange(2) > dataRange(2)
        endWaveform = round(templateLength - (spikeRange(2) - dataRange(2)));
        spikeRange(2) = round(dataRange(2));
      else
        endWaveform = templateLength;
      end
      iWaveform(startWaveform:endWaveform) = ones(1,numel(startWaveform:endWaveform));
      iWaveformFull = repmat(iWaveform,templateCh,1);
      iWaveformFull = reshape(iWaveformFull',1,templateLength,templateCh);
      if isnan(nTemplateWaveformsNew{filePointer}(iCluIDs,1,1))
        nTemplateWaveformsNew{filePointer}(iCluIDs,:,:) = iWaveformFull;
        nAmplitudes{filePointer}(iCluIDs) = 1;
      else
        nTemplateWaveformsNew{filePointer}(iCluIDs,:,:) = nTemplateWaveformsNew{filePointer}(iCluIDs,:,:) + iWaveformFull;
        nAmplitudes{filePointer}(iCluIDs) = nAmplitudes{filePointer}(iCluIDs) + 1;
      end
      waveform(:,startWaveform:endWaveform) = dat(logical(connected),(spikeRange(1):spikeRange(2))-(dataRange(1)-1));
      waveform = reshape(waveform',1,templateLength,templateCh);
      if isnan(templateWaveformsNew{filePointer}(iCluIDs,1,1))
        templateWaveformsNew{filePointer}(iCluIDs,:,:) = waveform;
        amplitudes{filePointer}(iCluIDs) = ampOI(i);
      else
        templateWaveformsNew{filePointer}(iCluIDs,:,:) = templateWaveformsNew{filePointer}(iCluIDs,:,:) + waveform;
        amplitudes{filePointer}(iCluIDs) = amplitudes{filePointer}(iCluIDs) + ampOI(i);
      end
    end
  else
    break
  end
  chunkInd = chunkInd+1;
end
for iFile = 1:nFiles
  waveforms{iFile} = templateWaveformsNew{iFile}./nTemplateWaveformsNew{iFile};
  amplitudes{iFile} = amplitudes{iFile}./nAmplitudes{iFile};
  if merge
    break
  end
end


%% Pick the largest waveforms
for iFile = 1:nFiles
  maxWaveforms{iFile} = zeros(numel(cluIDs), templateLength);
  maxChan{iFile} = zeros(size(cluIDs));
  chanMap{iFile} = zeros(numel(cluIDs), 3);
  if merge
    break
  end
end
for iFile = 1:nFiles
  for iWave = 1:numel(cluIDs)
%     prevMaxVal = 0;
    prevMaxDif = 0;
    for iChan = 1:templateCh
%       maxVal = max(abs(waveforms{iFile}(iWave,spikeCentreIndex-5:spikeCentreIndex+5,iChan)));
%       if maxVal > prevMaxVal
%         prevMaxVal = maxVal;
%         maxWaveforms{iFile}(iWave,:) = waveforms{iFile}(iWave,:,iChan);
%         maxChan{iFile}(iWave) = iChan;
% %         valChan = iChan;
%       end
      maxDif = max(abs(waveforms{iFile}(iWave,spikeCentreIndex-3:spikeCentreIndex+3,iChan) - mean(waveforms{iFile}(iWave,[1:spikeCentreIndex-10 spikeCentreIndex+10:end],iChan))));
      if maxDif > prevMaxDif
        prevMaxDif = maxDif;
        maxWaveforms{iFile}(iWave,:) = waveforms{iFile}(iWave,:,iChan);
        maxChan{iFile}(iWave) = iChan;
%         difChan = iChan;
      end
%       if iChan == templateCh
%         valChan
%         difChan
%       end
    end
    chanMap{iFile}(iWave,1) = cluIDs(iWave);
    chanMap{iFile}(iWave,2) = maxChan{iFile}(iWave);
    if cgs(cids == cluIDs(iWave)) == 0  % it's noise
      chanMap{iFile}(iWave, 3) = 0;
    elseif cgs(cids == cluIDs(iWave)) == 1 % it's MUA
      chanMap{iFile}(iWave, 3) = 1;
    else % it's a unit
      chanMap{iFile}(iWave, 3) = cluIDs(iWave);
    end
  end
  if merge
    break
  end
end


%% Clean-up waveforms
% for iFile = 1:1
%   waveExists = logical(chanMap{iFile}(:,2));
%   waveforms{iFile} = waveforms{iFile}(waveExists,:,:);
%   maxWaveforms{iFile} = maxWaveforms{iFile}(waveExists,:);
%   amplitudes{iFile} = amplitudes{iFile}(waveExists);
%   maxChan{iFile} = maxChan{iFile}(waveExists);
%   chanMap{iFile} = chanMap{iFile}(waveExists,:);
% end

for iFile = 1:nFiles
  emptyCount = ones(1,numel(cluIDs));
  for iClu = 1:numel(cluIDs)
    if nAmplitudes{iFile}(iClu) < 300
      emptyCount(iClu) = 0;
    end
  end
  emptyCount = logical(emptyCount);
  waveforms{iFile}(~emptyCount,:,:) = [];
  amplitudes{iFile}(~emptyCount) = [];
  maxWaveforms{iFile}(~emptyCount,:) = [];
  maxChan{iFile}(~emptyCount) = [];
  chanMap{iFile}(~emptyCount,:) = [];
end


%% Save waveforms
waveforms_temp = waveforms;
maxWaveforms_temp = maxWaveforms;
amplitudes_temp = amplitudes;
maxChan_temp = maxChan;
chanMap_temp = chanMap;
if merge
  iOutputFile = [dataDir filesep outputFile];
  datFile = datFileList;
  waveforms = waveforms_temp{1};
  maxWaveforms = maxWaveforms_temp{1};
  amplitudes = amplitudes_temp{1};
  maxChan = maxChan_temp{1};
  chanMap = chanMap_temp{1};
  save(iOutputFile, 'datFile','cluIDs','spikeCentreIndex','waveforms','maxWaveforms','amplitudes','maxChan','chanMap', '-v7.3');
else
  for iFile = 1:nFiles
    pathStr = fileparts(datFileList{iFile});
    iOutputFile = [pathStr filesep outputFile];
    datFile = datFileList{iFile}; %#ok<*NASGU>
    waveforms = waveforms_temp{iFile};
    maxWaveforms = maxWaveforms_temp{iFile};
    amplitudes = amplitudes_temp{iFile};
    maxChan = maxChan_temp{iFile};
    chanMap = chanMap_temp{iFile};
    save(iOutputFile, 'datFile','cluIDs','spikeCentreIndex','waveforms','maxWaveforms','amplitudes','maxChan','chanMap', '-v7.3');
  end
end
fclose all;


%% Display waveforms
if display
  displayWaveforms(cluIDs, maxWaveforms_temp, datFileList);
end