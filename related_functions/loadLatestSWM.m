  %% loading/compiling function
  
  load ~/GIT/SWM/compiled/latest
  revNumImp=importdata('/home/mrphys/bargip/GIT/SWM/version');
  revNumSC=regexp(revNumImp,'(r\S+-\d+)-','tokens');
  if isempty(revNumSC{1})
    revNumSC=regexp(revNumImp,'(r\S+)','tokens');
  end
  revNumSC=strrep(revNumSC{1}{1}{1},'-','_');
  if ~strcmp(revNumSC,revNum)
    revNum=revNumSC;
    temp=pwd;
    cd ~/GIT/SWM/compiled
    batchid=['bg_SWM_' revNum];
    compiledfun=qsubcompile('bg_SWM','batchid',batchid);
    compiledfun.revNum=revNum;
    fName=fullfile('~/GIT/SWM/compiled',batchid);
    save(fName,'compiledfun','revNum');
    save('~/GIT/SWM/compiled/latest','compiledfun','revNum')
    cd(temp)
  end
  disp(['loaded SWM function version ''' revNum ''''])