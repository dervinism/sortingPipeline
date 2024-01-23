% Input: directory where KiloSort output resides
% Output: clu and res in the legacy format (i.e. time is in samples), clu=0 is noise, clu=1 is MUA
function [clu, res] = getCluResFromKSdir(ksdir)

sp = loadKSdir(ksdir);

res = sp.st * sp.sample_rate;
u_clu = sort(unique(sp.clu), 'ascend');

if ~isempty(setdiff(u_clu, sp.cids))
  error(['sp.clu has units not in sp.cids: ' num2str(setdiff(u_clu, sp.cids))])
end
u_clu = sp.cids';

if u_clu(1) == 0 % convert template 0 to the next available id, as 0 will be reserved for noise clusters
  sp.clu(sp.clu == 0) = u_clu(end) + 1;
  u_clu = [u_clu(2:end); u_clu(end) + 1]; % new id
  sp.cgs = circshift(sp.cgs, -1); % move the classification of cluster 0 to the end, to correspond to u_clu
end

if u_clu(1) == 1 % convert template 1 to the next available id, as 0 will be reserved for noise clusters
  sp.clu(sp.clu == 1) = u_clu(end) + 1;
  u_clu = [u_clu(2:end); u_clu(end) + 1]; % new id
  sp.cgs = circshift(sp.cgs, -1); % move the classification of cluster 1 to the end, to correspond to u_clu
end

clu = sp.clu;

for i = 1:numel(u_clu)
  if sp.cgs(i) == 0 % noise cluster
    clu(clu == u_clu(i)) = 0;
  elseif sp.cgs(i) == 1 % MUA cluster
    clu(clu == u_clu(i)) = 1;
  end
end

clu = [numel(unique(clu)); clu]; % In the legacy format clu has the total number of clusters as its first element




