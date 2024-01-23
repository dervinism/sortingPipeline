% Performs postprocessing steps of the spike-sorting pipeline.
% Should be called after the manual refinement step in phy
% At this stage the postprocessing includes two steps:
% 1. computing the average spikewaveforms of the units (by averageing the raw recorded voltage traces)
% 2. calculating the quality of each unit
% Inputs: binaryFilename - the full name of the binary file (on which Kilosort ran)
%              noAIchans - number of AI (analog input) channels, which are not about the silicon probe (typically the last channel)
%                     sr - sampling rate of the binary file (in Hz, 3e4 by default)
function postprocessingPipeline(binaryFilename, noAIchans, sr)
if nargin < 3
  sr = 3e4;
end

[inp.dataDir, inp.dataFile, ext] = fileparts(binaryFilename);
inp.dataFile = [inp.dataFile, ext];
inp.display = false;
inp.merge   = true;
inp.chansIgnore = noAIchans;

extractWaveform(inp); % will extract the waveform of each unit, and save the results in waveforms.mat in the same directory
createQualityFileKilosort(inp.dataDir, sr); % This function computes autocorrelograms, and relies on CCG function (e.g. in R:\CSN\Shared\Dynamics\Code\github_cortex-lab_spikes\analysis\helpers\ )
