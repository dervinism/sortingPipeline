function [clu, res, templates, chanMap] = resCluFromKilosort(dirname, shankOI, shCh, chOI, probeFile)
% [clu, res, templates] = resCluFromKilosort(dirname, shankOI, shCh, chOI)
%
% This function extracts clu, res and templates variables from a kilosort
% output directory. It is a helper function to AnPSD_load,
% loadAsMUA_noResClu, loadAsRasterSparse, and extractWaveform.
%
% Inputs: dirname - a kilosort output directory containing npy files.
%         shankOI - the shank of interest.
%         shCh - a number of recording channels per shank.
%         chOI - channels of interest (so that you can look at specific
%                shank sections).
%         probeFile - a full path to a probe configuration file (forPRB*).
%                If you used kilosort2 to spikesort your data, you must
%                supply a waveform file. Otherwise, channel locations of
%                units will not be identified correctly.
%
% Outputs: clu - a vector with spike IDs. The first element of the vector
%                is the total number of unique units (MUA's are grouped
%                into a cluster ID 1).
%          res - spike index vector. Divide by the sampling frequency in
%                order to obtain the spike times. The following is true:
%                numel(res) == numel(clu)-1.
%          templates is templates = clu(2:end).
%          chanMap - unit channel locations. First column of the matrix
%                    contain cluster IDs, the second one contains probe
%                    channels, while the third one is the cluster type:
%                    0 - noise, 1 - mua, or a cluster ID if it is a single
%                    unit.

if nargin < 5
  probeFile = [];
end

ySelection = [];

sp = loadKSdir(dirname);

clu = sp.clu; %readNPY([dirname filesep 'spike_clusters.npy']);
res = sp.st * sp.sample_rate; %readNPY([dirname filesep   'spike_times.npy']); 
tmpl = sp.spikeTemplates+1;% readNPY([dirname filesep   'spike_templates.npy']); tmpl = tmpl+1; % they start from 0 (python way)
templateWaveforms = sp.temps; %readNPY([dirname filesep   'templates.npy']); % templates x time x channel
ycoordsCh = sp.ycoords;
xcoordsCh = sp.xcoords;

if isempty(clu) && isempty(res) && isempty(tmpl)
  templates = [];
  return
else
  assert(numel(clu) == numel(res) && numel(res) == numel(tmpl) && max(tmpl) <= size(templateWaveforms, 1))
end

cids = sp.cids; cgs = sp.cgs; %[cids, cgs] = readClusterGroupsCSV([dirname filesep 'cluster_groups.csv']);
%[cids, cgs] = readClusterGroupsCSV([dirname filesep 'cluster_group.tsv']);

uClu = double(unique(clu));

assert(max(abs(double(sort(torow(uClu))) - double(sort(torow(cids))))) == 0, 'should be fully compatible')
assert(~any(cgs >= 3), 'unsorted units remain')

if ~isempty(ySelection)
  [~, max_site] = max(max(abs(sp.temps),[],2),[],3); % the maximal site for each template
  spike_ycoord = sp.ycoords(max_site(sp.spikeTemplates+1));
  clu = clu(spike_ycoord >= ySelection(1) & spike_ycoord <= ySelection(2));
  res = res(spike_ycoord >= ySelection(1) & spike_ycoord <= ySelection(2));
  uClu = double(unique(clu));
end
clear max_site spike_ycoord

% Make sure no unit is named 0 or 1
if any(uClu == 0)
  m = max(uClu) + 1;
  clu(clu == 0) = m;
  cids(cids == 0) = m;
  uClu(uClu == 0) = m;
end
if any(uClu == 1)
  m = max(uClu) + 1;
  clu(clu == 1) = m;
  cids(cids == 1) = m;
  uClu(uClu == 1) = m;
end
templates = clu;
templateWaveforms2D = reshape(templateWaveforms, size(templateWaveforms, 1), []);
sh = zeros(size(clu)); % will hold the shank on which each spike resides (according to the template it's assigned to)
oi = zeros(size(clu)); % will be 1 if resides on the channel of interest

if ~isempty(probeFile)
  load(probeFile, 'ycoords','xcoords')
  ycoordsUnique = unique(ycoords);
  ycoordsCount = zeros(size(ycoords));
  for iCoord = 1:numel(ycoordsUnique)
    ycoordsCount(ycoords == ycoordsUnique(iCoord)) = sum(ycoords == ycoordsUnique(iCoord));
  end
  if numel(ycoords) <= 65
    ycoordsCount = ycoordsCount./min(ycoordsCount(ycoordsCount > 0));
  end
end

chanMap = [];
for u = torow(uClu)
  h = histc(tmpl(clu == u), 1:size(templateWaveforms, 1)); h = h/sum(h); %#ok<HISTC>
  w = squeeze(reshape(torow(h)*templateWaveforms2D, size(templateWaveforms, 2), size(templateWaveforms, 3)));
  % approximation for the average waveform of this unit. The computation is
  % done this way because some units can be assigned to more than one
  % template (e.g. after merge in phy)
  
  
  [~, pos] = max(abs(w(:))); if numel(pos) > 1; pos = pos(1); end
  % ceil(pos / size(tmp, 1)) -- the channel with the highest spike template waveform
  if isempty(probeFile)
    sh(clu == u) = ceil(ceil(pos / size(w, 1))/shCh);
    fprintf('Unit %d, on ch %d(%d) ==> it''s on shank %d and is ', u, ceil(pos / size(w, 1)), ceil(pos / size(w, 1))-1, ...
      ceil(ceil(pos / size(w, 1))/shCh))
    oi(clu == u) = sum(ceil(pos / size(w, 1)) == chOI);
    chanMap(end+1, :) = [u ceil(pos / size(w, 1)) NaN]; %#ok<AGROW>
    if cgs(cids == u) == 0     % it's noise
      clu(clu == u) = 0;
      chanMap(end, 3) = 0;
      fprintf('noise')
    elseif cgs(cids == u) == 1 % it's MUA
      clu(clu == u) = 1;
      chanMap(end, 3) = 1;
      fprintf('MUA')
    else
      chanMap(end,3) = u;
      fprintf('good')
    end
  else
    pos = ceil(pos / size(w, 1));
    ycoordCh = ycoordsCh(pos);
    if exist('ycoordsCount','var')
      ycoordChCount = ycoordsCount(pos);
    else
      ycoordChCount = 1;
    end
    xcoordCh = xcoordsCh(pos);
    posY = find(ycoordCh == ycoords);
    posX = find(xcoordCh == xcoords(posY));
    pos = posY(1) + floor((posX - 1)/ycoordChCount)*shCh;
    sh(clu == u) = ceil(pos/shCh);
    fprintf('Unit %d, on ch %d(%d) ==> it''s on shank %d and is ', u, pos, pos-1, ceil(pos/shCh))
    oi(clu == u) = sum(pos == chOI);
    chanMap(end+1, :) = [u pos NaN]; %#ok<AGROW>
    if cgs(cids == u) == 0     % it's noise
      clu(clu == u) = 0;
      chanMap(end, 3) = 0;
      fprintf('noise')
    elseif cgs(cids == u) == 1 % it's MUA
      clu(clu == u) = 1;
      chanMap(end, 3) = 1;
      fprintf('MUA')
    else
      chanMap(end,3) = u;
      fprintf('good')
    end
  end
  fprintf('\n')
end

if shankOI
  clu = [numel(unique(clu(sh == shankOI & oi))); clu(sh == shankOI & oi)];
  res = res(sh == shankOI & oi);
  templates = templates(sh == shankOI & oi);
else
  clu = [numel(unique(clu(oi))); clu(oi)];
  res = res(oi);
  templates = templates(oi);
end
clu = round(clu);
res = round(res);
templates = round(templates);