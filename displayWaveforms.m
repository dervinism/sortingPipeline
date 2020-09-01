function displayWaveforms(cluIDs, maxWaveforms, datFileList)

if (iscell(datFileList) && numel(datFileList) > 1) && (iscell(maxWaveforms) && numel(maxWaveforms) > 1)
  nFiles = numel(datFileList);
else
  nFiles = 1;
  if iscell(datFileList) && numel(datFileList) > 1
    clear datFileList
    datFileList{1} = 'all';
  elseif ~iscell(datFileList)
    datFileList_temp = datFileList;
    clear datFileList
    datFileList{1} = datFileList_temp;
  end
  if ~iscell(maxWaveforms)
    maxWaveforms_temp = maxWaveforms;
    clear maxWaveforms
    maxWaveforms{1} = maxWaveforms_temp;
  end
end

for iFile = 1:nFiles
  for iWave = 1:numel(cluIDs) %#ok<*UNRCH>
    figure('units', 'normalized', 'position', [0.002, .04, 1, .88], 'Visible','on');
    plot(maxWaveforms{iFile}(iWave,:))
    title(['file ' datFileList{iFile} '   unit ' num2str(cluIDs(iWave))], 'Interpreter','none');
    ylabel('\muV');
    xlabel('data points')
  end
end