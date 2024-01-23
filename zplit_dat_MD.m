% Creates .clu, .res, and .tmpl files corresponding to each individual .dat
% of those ztitched into the compound .dat file
% Input:
%     datFilename - name of the compound .dat file (must end with '.dat'),
%                   the file itself need not exist (it's the corresponding .mat file that is really needed)
%     dropFet     - (optional, false by default) If true will not do anything with fet files  

function zplit_dat_MD(datFilename, dropFet)
if nargin < 1 || ~strcmp(datFilename(end-3:end), '.dat')
  error('wrong inputs');
end

if ~exist([datFilename(1:end-4) '.mat'], 'file')  
  error([datFilename(1:end-4) '.mat not exists?!']);
end

if nargin < 2
  dropFet = false;
end

ztitchInfo = load([datFilename(1:end-4) '.mat']);
assert(length(ztitchInfo.datFileList) == length(ztitchInfo.dataPoints));

for sh = 1:100000 % Loop over shanks
  if ~exist([datFilename(1:end-4) '.clu.' num2str(sh)], 'file')  
    break; % we're done
  end
  disp(['Shank ' num2str(sh) ' ----- ']);
  if ~exist([datFilename(1:end-4) '.res.' num2str(sh)], 'file')  
    dlmwrite([datFilename(1:end-4) '.res.' num2str(sh)], ...
      hdf5read([datFilename(1:end-4) '.kwik'], ...
      ['/channel_groups/' num2str(sh) '/spikes/time_samples']), ...
      'precision', '%15.0f');
  end
  
  clu = load([datFilename(1:end-4) '.clu.' num2str(sh)]);
  res = load([datFilename(1:end-4) '.res.' num2str(sh)]);
  templates = load([datFilename(1:end-4) '.tmpl.' num2str(sh)]);
  assert(numel(clu) == numel(res)+1);
  clu = clu(2:end);
  
  [fileDir, fileName] = fileparts(datFilename(1:end-4));
  if ~dropFet && exist([datFilename(1:end-4) '.fet.' num2str(sh)], 'file')
    fH = fopen([datFilename(1:end-4) '.fet.' num2str(sh)], 'r');
  elseif ~dropFet && exist([fileDir filesep '_klustakwik' filesep fileName '.fet.' num2str(sh)], 'file')
    fH = fopen([fileDir filesep '_klustakwik' filesep fileName '.fet.' num2str(sh)]);
  else
    fH = -1;
  end
  clear fileDir fileName;
    
  if fH > 0
    dim = textscan(fH, '%d', 1); dim = dim{1};
    fet = textscan(fH, '%d', dim*length(res));
    fclose(fH);
    fet = reshape(fet{1}, dim, []);
    assert(size(fet, 2) == length(res));  
  else
    if ~dropFet
      disp('No fet file...');
    end
    fet = [];
  end
  
  for z = 1:length(ztitchInfo.dataPoints) % loop over the original .dat files
    disp([ztitchInfo.datFileList{z} ' ...']);
    I = res <= ztitchInfo.dataPoints(z);
    newClu = clu(I); newClu = [length(unique(newClu)); newClu]; %#ok<*AGROW>
    newRes = res(I);
    newTemplates = templates(I);
    dlmwrite([ztitchInfo.datFileList{z}(1:end-4) '.clu.' num2str(sh)], ...
      newClu);
    dlmwrite([ztitchInfo.datFileList{z}(1:end-4) '.res.' num2str(sh)], ...
      newRes, 'precision', '%15.0f');
    dlmwrite([ztitchInfo.datFileList{z}(1:end-4) '.tmpl.' num2str(sh)], ...
      newTemplates);
    clu = clu(~I);
    res = res(~I) - ztitchInfo.dataPoints(z);
    templates = templates(~I);
    
    if ~isempty(fet)
      newFet = fet(:,I);
      dlmwrite([ztitchInfo.datFileList{z}(1:end-4) '.fet.' num2str(sh)], dim);
      dlmwrite([ztitchInfo.datFileList{z}(1:end-4) '.fet.' num2str(sh)], newFet', '-append', 'delimiter', ' ');
      fet = fet(:, ~I);
    end
    
  end
end
