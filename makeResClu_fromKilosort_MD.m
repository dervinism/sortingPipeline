cd(fileparts(mfilename('fullpath'))) % change to the directory where this script resides

dirname = 'R:\CSN\Shared\Dynamics\Data\M190523_B_MD\201906281';
datFilename = 'stitchedData'; % no .dat in the end needed!
ySelection = []; % if we want to use only a portion of the recording, restricted to a certain part of the probe
nCh = 15;
chsh = 15;

if ispc % for readNPY, readClusterGroupsCSV
  addpath(genpath('..\..\github_kwikteam_npy-matlab'))
  addpath(genpath('..\..\github_cortex-lab_spikes'))
else % linux machine
  addpath(genpath('/data/nick/code/npy-matlab'))
  addpath(genpath('/data/mush/github_cortex-lab_spikes'))
end
% addpath('/data/mush/mize') % for computeMUSW

if ~exist('ySelection', 'var')
  ySelection = [];
end


sp = loadKSdir(dirname); % (by default) this excludes spikes designated as 'noise' 

clu = sp.clu; %readNPY([dirname filesep 'spike_clusters.npy']);
res = sp.st * sp.sample_rate; %readNPY([dirname filesep   'spike_times.npy']); 
tmpl = sp.spikeTemplates+1;% readNPY([dirname filesep   'spike_templates.npy']); tmpl = tmpl+1; % they start from 0 (python way)
templateWaveforms = sp.temps; %readNPY([dirname filesep   'templates.npy']); % templates x time x channel

assert(numel(clu) == numel(res) && numel(res) == numel(tmpl) && max(tmpl) <= size(templateWaveforms, 1))

cids = double(sp.cids); cgs = sp.cgs; %[cids, cgs] = readClusterGroupsCSV([dirname filesep 'cluster_groups.csv']);
%[cids, cgs] = readClusterGroupsCSV([dirname filesep 'cluster_group.tsv']);

uClu = double(unique(clu));

assert(max(abs(sort(torow(uClu)) - sort(torow(cids)))) == 0, 'should be fully compatible')
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
templateWaveforms2D = reshape(templateWaveforms, size(templateWaveforms, 1), []);
sh = zeros(size(clu)); % will hold the shank on which each spike resides (according to the template it's assigned to)
templates = clu;

chanMap = [];
for u = torow(uClu)
  h = histc(tmpl(clu == u), 1:size(templateWaveforms, 1)); h = h/sum(h);
  w = squeeze(reshape(torow(h)*templateWaveforms2D, size(templateWaveforms, 2), size(templateWaveforms, 3)));
  % approximation for the average waveform of this unit. The computation is
  % done this way because some units can be assigned to more than one
  % template (e.g. after merge in phy)
  
  
  [~, pos] = max(abs(w(:))); if numel(pos) > 1; pos = pos(1); end
  % ceil(pos / size(tmp, 1)) -- the channel with the highest spike template waveform
  sh(clu == u) = ceil(ceil(pos / size(w, 1))/chsh);
  fprintf('Unit %d, on ch %d(%d) ==> it''s on shank %d and is ', u, ceil(pos / size(w, 1)), ceil(pos / size(w, 1))-1, ...
    ceil(ceil(pos / size(w, 1))/chsh))
  chanMap(end+1, :) = [u ceil(pos / size(w, 1)) NaN];
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
  fprintf('\n')
end

save([dirname filesep datFilename '.chanMap.mat'], 'chanMap') % on which channel each unit resides (i.e. has largest spikes), and what type it is 0[noise]/1[mua]/[good]
% chanMap has exactly the same units as spike_clusters.npy, except for units 0,1 which have different meaning and will be renamed
for i = 1:nCh/chsh
  dlmwrite([dirname filesep datFilename '.res.' num2str(i)], res(sh == i), 'precision', '%15.0f');
  dlmwrite([dirname filesep datFilename '.clu.' num2str(i)], [numel(unique(clu(sh == i))); clu(sh == i)], 'precision', '%15.0f');
  dlmwrite([dirname filesep datFilename '.tmpl.' num2str(i)], templates(sh == i), 'precision', '%15.0f');
end
return

dlmwrite([dirname filesep datFilename '.res.99' ], res, 'precision', '%15.0f');
dlmwrite([dirname filesep datFilename '.clu.99' ], [numel(unique(clu)); clu], 'precision', '%15.0f');
dlmwrite([dirname filesep datFilename '.tmpl.99' ], templates, 'precision', '%15.0f');
[avW, clN] = computeMUSW([dirname filesep datFilename], 99, nCh);
delete([dirname filesep datFilename '.res.99']) % we did that to compute all the waveforms in one go...
delete([dirname filesep datFilename '.clu.99'])
delete([dirname filesep datFilename '.tmpl.99'])

for i = 1:nCh/chsh
  %unique(clu(sh == i))
  avWaveform = avW(ismember(clN, unique(clu(sh == i))), :, :);
  clusterNumbers = clN(ismember(clN, unique(clu(sh == i))));
  save([dirname filesep datFilename '.musw.' num2str(i) '.mat'], 'avWaveform', 'clusterNumbers')
end
