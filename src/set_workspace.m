%   Project Title: SET WORKSPACE FOR SCRIPTS
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220810.0
%   Previous Version: n/a
%   Summary: this script is an initializer and workspace variable setup for
%   all scripts in this repository

%% SET PARAMETERS
%## TIME
TT = tic;

%## folders
if ~exist('ADD_ALL_SUBMODS','var')
    ADD_ALL_SUBMODS = false;
end

%##
SUBMODS_DIR = [fileparts(SRC_DIR) filesep 'submods'];
FUNCS_DIR = [fileparts(fileparts(SRC_DIR)) filesep 'MindInMotion_Functions' filesep 'src'];

%##
path(path,SUBMODS_DIR)
path(path,SRC_DIR)
path(path,FUNCS_DIR);
%% HARDCODE PATHS STRUCT =============================================== %%
PATHS = [];
RMV_MODS = {'dipfit'};
if ADD_ALL_SUBMODS
    SUBMODULES_GENPATH = {'cleanline'};
    tmpd = dir(SUBMODS_DIR);
    tmpd = tmpd(3:end);
    SUBMODULES = {tmpd.name};    
else
    SUBMODULES = {'fieldtrip','eeglab','sift','postamicautility',...
            'iclabel','viewprops','powpowcat', ...
            'bemobil_pipeline','gait_tracking_w_imu','eeglab_specparam'};
    SUBMODULES_GENPATH = {};
end
%--
for i = 1:length(RMV_MODS)
    ind = strcmp(RMV_MODS{i},SUBMODULES);
    if any(ind)
        SUBMODULES{ind} = [];
    end
end
SUBMODULES = SUBMODULES(~cellfun(@isempty,SUBMODULES));
%(03/27/2025) JS, removing dipfit as it leads to odd behavior 'dipfit'


%## add submodules
if ispc
    DELIM = ';';
else
    DELIM = ':';
end
PATHS.PATHS = cell(length(SUBMODULES),1);
for ss = 1:length(SUBMODULES)
    if any(strcmp(SUBMODULES{ss},SUBMODULES_GENPATH))
        fprintf('Adding submodule using genpath(): %s...\n',[SUBMODS_DIR filesep SUBMODULES{ss}]);
        a_ftmp = unix_genpath([SUBMODS_DIR filesep SUBMODULES{ss}]);
        a_ftmp = split(a_ftmp,DELIM); a_ftmp = a_ftmp(~cellfun(@isempty,a_ftmp));
        cellfun(@(x) path(x,path),a_ftmp);
        cellfun(@(x) fprintf('Adding functions in: %s...\n',x),a_ftmp);
    else
        fprintf('Adding submodule: %s...\n',[SUBMODS_DIR filesep SUBMODULES{ss}]);
        path([SUBMODS_DIR filesep SUBMODULES{ss}],path);
    end
    %-- spm12 exception
    if strcmp(SUBMODULES{ss},'spm12') && exist([SUBMODS_DIR filesep SUBMODULES{ss} filesep SUBMODULES{ss}],'dir')
        path(path,[SUBMODS_DIR filesep SUBMODULES{ss} filesep SUBMODULES{ss}]);
    end
    %(02/26/2025) JS, this could defintely use a robustness improvement
    PATHS.PATHS{ss} = [SUBMODS_DIR filesep SUBMODULES{ss}];
end
%## special paths
% %- current dir
% PATHS.curr_dir = TMP_PWD;
%- submods path
PATHS.submods_dir = SUBMODS_DIR;
%- src folder
PATHS.src_dir = SRC_DIR;
%- _data folder
PATHS.data_dir = fullfile(fileparts(fileparts(SRC_DIR)),'_data');
%- _functions folder
PATHS.functions_dir = FUNCS_DIR;
a_ftmp = unix_genpath(PATHS.functions_dir);
a_ftmp = split(a_ftmp,DELIM);
a_ftmp = a_ftmp(~cellfun(@isempty,a_ftmp));
cellfun(@(x) path(x,path),a_ftmp);
cellfun(@(x) fprintf('Adding functions in: %s...\n',x),a_ftmp);
a_ftmp = char.empty;
%(04/27/2025) JS, making sure that paths are added to the top of the path
%by using the form path(path_to_add,path);
% ----------------------------------------------------------------------- %
%% ADDPATH for FIELDTRIP Toolboxbemobil
if contains('fieldtrip',SUBMODULES,'IgnoreCase',true)
    fprintf('Starting fieldtrip...\n');
    ft_defaults;
end
%% INITIALIZE MIM & EEGLAB
%start EEGLAB if necessary
if contains('sift',SUBMODULES,'IgnoreCase',true) %&& ~contains('cleanline',SUBMODULES,'IgnoreCase',true)
    fprintf('Starting SIFT...\n');
    StartSIFT;
end
%- always start eeglab last.
ALLCOM=cell.empty;
ALLEEG=struct.empty;
STUDY=struct.empty;
CURRENTSET=double.empty;
CURRENTSTUDY=double.empty;
eeglab;
if contains('eeglab',SUBMODULES,'IgnoreCase',true)
    fprintf('Starting EEGLAB...\n');
    %- EEGLAB folder
    PATHS.eeglab_dir = [SUBMODS_DIR filesep 'eeglab'];
    eeglab nogui;
    %(02/26/2025) JS, adding "nogui" keyword to see if improves scripting
    %performance
end

%% PARPOOL SETUP ======================================================= %%
if ~ispc
    pop_editoptions('option_storedisk',1, ...
        'option_savetwofiles',1, ...
        'option_single',1, ...
        'option_memmapdata',0, ...
        'option_computeica',0, ...
        'option_saveversion6',1, ...
        'option_scaleicarms',1, ...
        'option_rememberfolder',1);
    disp(['SLURM_JOB_ID: ', getenv('SLURM_JOB_ID')]);
    disp(['SLURM_CPUS_ON_NODE: ', getenv('SLURM_CPUS_ON_NODE')]);
    %## allocate slurm resources to parpool in matlab
    %- get cpu's on node and remove a few for parent script.
    SLURM_POOL_SIZE = str2double(getenv('SLURM_CPUS_ON_NODE'));
    %- create cluster
    pp = parcluster('local');
    %- Number of workers for processing (NOTE: this number should be higher
    %then the number of iterations in your for loop)
    fprintf('Number of workers: %i\n',pp.NumWorkers);
    fprintf('Number of threads: %i\n',pp.NumThreads);
    %- make meta data dire1ory for slurm
    if exist(SRC_DIR,'var')
        mkdir([SRC_DIR filesep '_slurm_scratch' filesep getenv('SLURM_JOB_ID')])
        pp.JobStorageLocation = [SRC_DIR filesep '_slurm_scratch' filesep getenv('SLURM_JOB_ID')];
    else
        fprintf('Slurm scratch location is not defined...\n\n');
    end
    %- create your p-pool (NOTE: gross!!)
    pPool = parpool(pp, SLURM_POOL_SIZE, 'IdleTimeout', Inf);
else
    pop_editoptions('option_storedisk',1, ...
        'option_savetwofiles',1, ...
        'option_single',1, ...
        'option_memmapdata',0, ...
        'option_computeica',0, ...
        'option_saveversion6',1, ...
        'option_scaleicarms',1, ...
        'option_rememberfolder',1);
    SLURM_POOL_SIZE = 1;
end
%% FINAL PRINTS
fprintf(1,'Using pathing:\n-WORKSPACE: %s\n-SUBMODULES: %s\n-FUNCTIONS: %s\n-DATA: %s\n',PATHS.src_dir,PATHS.submods_dir,PATHS.functions_dir,PATHS.data_dir);
%## TIME
fprintf('Done with workplace setup: %0.2f\n\n',toc(TT));
%% (SUBFUNCTIONS) ====================================================== %%
function [p] = unix_genpath(d)
    %GENPATH Generate recursive toolbox path.
    %   P = GENPATH returns a character vector containing a path name 
    %   that includes all the folders and subfolders below MATLABROOT/toolbox, 
    %   including empty subfolders.
    %
    %   P = GENPATH(FOLDERNAME) returns a character vector containing a path 
    %   name that includes FOLDERNAME and all subfolders of FOLDERNAME, 
    %   including empty subfolders.
    %   
    %   NOTE 1: GENPATH will not exactly recreate the original MATLAB path.
    %
    %   NOTE 2: GENPATH only includes subfolders allowed on the MATLAB
    %   path.
    %
    %   See also PATH, ADDPATH, RMPATH, SAVEPATH.

    %   Copyright 1984-2018 The MathWorks, Inc.
    %------------------------------------------------------------------------------

    % String Adoption
    if nargin > 0
        d = convertStringsToChars(d);
    end

    if nargin==0
      p = genpath(fullfile(matlabroot,'toolbox'));
      if length(p) > 1, p(end) = []; end % Remove trailing pathsep
      return
    end

    % initialise variables
    classsep = '@';  % qualifier for overloaded class directories
    packagesep = '+';  % qualifier for overloaded package directories
    p = '';           % path to be returned

    % Generate path based on given root directory
    files = dir(d);
    if isempty(files)
      return
    end

    % Add d to the path even if it is empty.
    p = [p d pathsep];

    % set logical vector for subdirectory entries in d
    isdir = logical(cat(1,files.isdir));
    %
    % Recursively descend through directories which are neither
    % private nor "class" directories.
    %
    dirs = files(isdir); % select only directory entries from the current listing

    for i=1:length(dirs)
       dirname = dirs(i).name;
       if    ~strcmp( dirname,'.')          && ...
             ~strcmp( dirname,'..')         && ...
             ~strncmp( dirname,classsep,1) && ...
             ~strncmp( dirname,packagesep,1) && ...
             ~strcmp( dirname,'private')    && ...
             ~strcmp( dirname,'resources') && ...
             ~strcmp( dirname,'__archive') && ...
             ~strcmp( dirname,'_compiles') && ...
             ~strcmp( dirname,'_mcc_funcs')
          p = [p genpath([d filesep dirname])]; % recursive calling of this function.
       end
    end
end
% %% ===================================================================== %%
% function dpath = finddir(start_dir,keyword)
%     %## initialise variables
%     classsep = '@';  % qualifier for overloaded class directories
%     packagesep = '+';  % qualifier for overloaded package directories
%     dpath = '';
% 
%     %## GET DIRNAMES
%     tmpd = dir(start_dir);
%     dirnames = {tmpd.name};
%     inds = ~strcmp( dirnames,'.') & ...
%         ~strcmp( dirnames,'..') & ...
%         ~strncmp( dirnames,classsep,1) & ...
%         ~strncmp( dirnames,packagesep,1) & ...
%         ~strcmp( dirnames,'private') & ...
%         ~strcmp( dirnames,'resources') & ...
%         ~strcmp( dirnames,'__archive') & ...
%         ~strcmp( dirnames,'_compiles');
%     tmpd = tmpd(inds);
%     fpaths = cellfun(@(x) [tmpd.folder filesep x],{tmpd.name},'UniformOutput',false); 
%     chk = strcmp({tmpd.name},keyword);
%     if any(chk)
%         dpath = [tmpd(chk).folder filesep tmpd(chk).name];
%     end
% end
% 
% %##
% function dpath = hlp_finddir(start_dir)
%     %## initialise variables
%     dpath = [];
%     if isempty(start_dir)
%         %## GET DIRNAMES
%         tmpd = dir(start_dir);
%         dirnames = {tmpd.name};
%         inds = ~strcmp( dirnames,'.') & ...
%             ~strcmp( dirnames,'..') & ...
%             ~strncmp( dirnames,classsep,1) & ...
%             ~strncmp( dirnames,packagesep,1) & ...
%             ~strcmp( dirnames,'private') & ...
%             ~strcmp( dirnames,'resources') & ...
%             ~strcmp( dirnames,'__archive') & ...
%             ~strcmp( dirnames,'_compiles');
%     end
%     tmpd = tmpd(inds);
%     fpaths = cellfun(@(x) [tmpd.folder filesep x],{tmpd.name},'UniformOutput',false); 
%     chk = strcmp({tmpd.name},keyword);
% end