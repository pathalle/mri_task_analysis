SPM.xX.X
figure;plot(SPM.xX.X(:,1))
figure;plot(SPM.xX.X(:,1:3))
figure;plot(SPM.xX.X(:,1:3));xlim([0 100])
figure;plot(SPM.xX.X(:,1:2));xlim([0 100])
corr(SPM.xX.X(:,1),SPM.xX.X(:,2))
corr(SPM.xX.X(:,1),SPM.xX.X(:,3))
corr(SPM.xX.X(:,1),SPM.xX.X(:,16))
corr(SPM.xX.X(:,1),SPM.xX.X(:,12))
[r,p]=corr(SPM.xX.X(:,1),SPM.xX.X(:,12))
SPM.xX.name
SPM.xX.name{1}
SPM.xX.Bcov
SPM.Sess
SPM.Sess(1).U
SPM.Sess(1).U.ons
SPM.Sess(1).U.u
figure;plot(SPM.xX.X(:,1:3))
figure;plot(SPM.xX.X(:,1:2))
figure;plot(SPM.xX.X(:,3:4))
subject = 'AllRead_pilot06';
logfile = cellstr(ls([paths.study, paths.logs, task '\' subject, '\*.txt']));
EVall = nan(128,1);
PEall = nan(128,1);
pGainAll = nan(128,1);
taskList = { 'learning_4_gfg' };
paths.study = 'O:\studies\allread\';
paths.analysis = 'analysis\';
paths.logs = 'logs\mri\Learning_Task\';
paths.pps = 'data\mri\Learning_Task\preprocessing\';
subjects = {};
for i=[ 7 ]
sub = sprintf('AllRead_pilot%02d',i);
subjects{end + 1} = sub;
end
% create batchfile for each subject
for i=1:length(subjects)
for t = 1:length(taskList)
% batches{length(taskList)*i+t-1} = create_glm_learning_4gfg(paths,taskList{t},subjects{i},0,3);
batches{length(taskList)*i+t-1} = create_glm_learning_4gfg(paths,taskList{t},subjects{i},0,3);
end
end
taskList = { 'learning_4_gfg' };
paths.study = 'O:\studies\allread\Piloting\';
paths.analysis = 'analysis\';
paths.logs = 'logs\mri\Learning_Task\';
paths.pps = 'data\mri\Learning_Task\preprocessing\';
% create batchfile for each subject
for i=1:length(subjects)
for t = 1:length(taskList)
% batches{length(taskList)*i+t-1} = create_glm_learning_4gfg(paths,taskList{t},subjects{i},0,3);
batches{length(taskList)*i+t-1} = create_glm_learning_4gfg(paths,taskList{t},subjects{i},0,3);
end
end
spm_jobman('interactive',batches{1})