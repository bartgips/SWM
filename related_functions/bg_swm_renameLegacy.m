function [ cfg ] = bg_swm_renameLegacy( cfg )
% [ cfg ] = bg_swm_renameLegacy( cfg )
% rename old structure arrays to new format

[~,cfg]=isfieldi(cfg,{'fName','varName'});

try
  cfg.numWindows=cfg.numTemplates;
  cfg=rmfield(cfg,'numTemplates');
end
try
  cfg.winPerTrial=cfg.numWin;
  cfg=rmfield(cfg,'numWin');
end

if isfield(cfg,'clust')
  for n=1:numel(cfg.clust)
    clust=cfg.clust{n};
    try
      clust.linIdx=clust.lidx;
      clust=rmfield(clust,'lidx');
    end
    try
      clust.tIdx=clust.tidx;
      clust=rmfield(clust,'tidx');
    end
    try
      clust.numWindows=clust.numTemplates;
      clust=rmfield(clust,'lidx');
    end
    cfg.clust{n}=clust;
  end
end

if isfield(cfg,'best_clust')
  if ~iscell(cfg.best_clust)
    clust=cfg.best_clust;
    cfg.best_clust={clust};
  end
  for n=1:numel(cfg.best_clust)
    clust=cfg.best_clust{n};
    try
      clust.linIdx=clust.lidx;
      clust=rmfield(clust,'lidx');
    end
    try
      clust.tIdx=clust.tidx;
      clust=rmfield(clust,'tidx');
    end
    try
      clust.numWindows=clust.numTemplates;
      clust=rmfield(clust,'lidx');
    end
    cfg.best_clust{n}=clust;
  end
end

