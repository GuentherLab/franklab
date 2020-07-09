
function varargout=fl_process(STEPS,varargin)

% FL_PROCESS Internal function : Frank Lab fMRI preprocessing/analysis procedures
%
%
% fl_process(STEPS,...) runs preprocessing/analysis steps
%     STEPS      : one or several among 'ALL', 'PREPROC', 'PREPROC.IMPORT','PREPROC.APPEND','MODEL', 'QA.CREATE', 'QA.PLOT'
% 
% **** PREPROC (or ALL) ****
% fl_process('PREPROC',...) imports functional/structural data and runs preprocessing steps
%     [,'subject_id'],id              : subject IDs = folder names where data will be stored
%     [,'subject_info',files]         : subject data information files (one .cfg file per subject; e.g. {'/data/subject1.cfg','/data/subject2.cfg',...})
%     [,'subject_info',filter]        : alternatively, path to subject data information files (e.g. '/data/subject*.cfg')
%     [,'subject_info',files]         : alternatively, path to project files to be copied (one .mat file per subject)
%     [,'pipeline_info',file]         : path to additional preprocessing pipeline information file, if any (e.g. 'example_preprocfile.cfg')
%     [,'parallel',flag]              : [1/0] set to true to run using cluster parallelization resources, or set to false to run locally (default: 0)
%     [,'overwrite',flag]             : [1/0] set to true to overwrite existing subject-id folders, or set to false to skip already-existing subject-id folders
%     [,'root_folder',folder]         : output folder (where subject-ID folder names will be created/stored; if not specified root_folder will be the default FL data folder; this field is disregarded if using explicit full paths as subject_id's)
%     [,'donotexpand',flag]           : [1/0] when entering a single subject_id string and multipel subject_info files, donotexpand=0 expands SUBJECT_ID to SUBJECT_ID/subject_info{1} SUBJECT_ID/subject_info{2} ... (a different project per subject) while donotexpand=1 expands SUBJECT_ID to SUBJECT_ID,1 SUBJECT_ID,2 ... (single project with multiple subjects) (default:0)
%
% example syntax:
%   fl_process PREPROC example subject_info /project/busplab/software/fl/example_subject.cfg pipeline_info /project/busplab/software/fl/example_preprocfile.cfg 
%
%
% note: alternative use 'PREPROC.IMPORT' (and do not enter a 'pipeline_info' field) to only import the data (without preprocessing it yet) 
% note: alternative use 'PREPROC.APPEND' (and do not enter a 'subject_info' field) to run additional preprocessing steps to an already preprocessed dataset
%
%
% **** MODEL (or ALL) ****
% fl_process('MODEL',...) runs first-level analysis steps
%     [,'subject_id'],id              : subject IDs = folder names where data is being stored
%     [,'design_info',file]           : path to additional first-level design information file, if any (e.g. 'example_modelfile.cfg')
%     [,'subject_info',files]         : (optional) alternative subject data information files containing subject-specific design information (one .cfg file per subject; e.g. {'/data/subject1.cfg','/data/subject2.cfg',...})
%
% example syntax:
%   fl_process MODEL example design_info /project/busplab/software/fl/example_modelfile.cfg 
%
%
% **** QA.CREATE (or ALL) ****
% fl_process('QA.CREATE',...) creates Quality Assurance plots
%     [,'subject_id'],id              : subject IDs = folder names where data is being stored
%     [,'qa_plist',plots]             : list of QA plots to create (see "help conn_qaplots" for additional details)
%     [,'qa_set',label]               : functional dataset (label or index) for preprocessing plots [1]
%
% example syntax:
%   fl_process QA.CREATE example
%
%
% **** QA.PLOTS ****
% fl_process('QA.PLOTS',...) displays Quality Assurance plots
%     [,'subject_id],id              : subject IDs = folder names where data is being stored
%
% example syntax:
%   fl_process QA.PLOTS example
%
%
% **** UTILS ****
% fl_process('UTILS_ROOTFOLDER') returns root data folder
% fl_process('UTILS_SETPERMISSIONS',folder) sets write permissions to group for all files within "folder"
%
% fl_process('OPEN',id) displays id dataset in CONN's gui
% fl_process('PARALLEL.REPORT',id) checks the status of submitted job
% fl_process('PARALLEL.DELETE',id) cancels/deletes submitted job
% fl_process('MODEL.PLOTS',subject_id) displays first-level contrasts
%
%
% notes on 'subject_id' field:
%   a)  if subject_id contains no path (or only relative path) information, subject_id will point to '/$ROOT_FOLDER/$SUBJECT_ID' (with $ROOT_FOLDER by default $FL_PATH/REPOSITORY)
%              e.g. fl_process PREPROC sub-0001 is equivalent to fl_process PREPROC /project/busplab/software/fl/REPOSITORY/sub-0001
%              e.g. fl_process PREPROC Pert/sub-0001 is equivalent to fl_process PREPROC /project/busplab/software/fl/REPOSITORY/Pert/sub-0001
%       if subject_id contains full path information, subject_id will point to the specified location
%              e.g. fl_process PREPROC /project/busplab/Pert/sub-0001 will import data in said directory/folder
%   b)  if not specified, subject ID will be derived from subject_info filenames 
%         e.g. fl_process PREPROC subject_info /data/DataSubject*.cfg  (matching several files DataSubject1.cfg DataSubject2.cfg ...) 
%              will use subject_id={'DataSubject1','DataSubject2',...};
%   c)  subject_id can also indicate a "format" string containing the pattern * or '%s', with the string being imported from the subject_info filenames
%         e.g. fl_process PREPROC subject_id 'Exp/*' subject_info /data/DataSubject*.cfg
%              will use subject_id={'Exp/DataSubject1','Exp/DataSubject2',...};
%   d)  subject_id can also indicate a "format" string containing the pattern '%d' (or '%0Nd'), with only the numeric part being now imported from the subject_info filenames
%         e.g. fl_process PREPROC subject_id 'Exp/sub-%04d' subject_info /data/DataSubject*.cfg
%              will use subject_id={'Exp/sub-0001','Exp/sub-0002',...};
%   e)  in datasets containing multiple subjects (e.g. data imported from CONN), use the format "mydataset,N" in the subject_id field to indicate the "N"th subject within the "mydataset" dataset
%         e.g. fl_process MODEL 'Exp/conn_project,%d' subject_info /data/DataSubject*.cfg design_info modelfile.cfg
%              will run first-level analyses for the subjects defined in the Exp/conn_project.mat CONN project file
%
%


persistent done_init rootfolder;

varargout=cell(1,max(1,nargout));
fileout=[];

if (ischar(STEPS)||(iscell(STEPS)&&numel(STEPS)==1))&&~ismember(lower(STEPS),{'all','preproc','preproc.append','preproc.import','model','qa.create','qacreate','qa.plot','qaplot','qa'})
    STEPS=char(STEPS);
    switch(lower(STEPS))
        case {'init','initforce'}
            if strcmpi(STEPS,'initforce') || isempty(done_init) || ~(isequal(done_init,0)||isequal(done_init,path))
                if isempty(which('evlab17'))
                    folderOptions_evlab={fullfile(fileparts(fileparts(which(mfilename))),'evlab17'),'/software/evlab17'};      % folders where to search for EVLAB17
                    tfilename=which('evlab17');
                    for n=1:numel(folderOptions_evlab),
                        if ~isempty(dir(folderOptions_evlab{n})),
                            if ~isempty(tfilename),
                                fprintf('warning: temporarily removing path entry %s\n  please remove this entry from your path to avoid this message in the future\n',fileparts(tfilename));
                                rmpath(fileparts(tfilename));
                            end
                            addpath(folderOptions_evlab{n});
                            break;
                        end
                    end
                    if isempty(which('evlab17')), error('please add path to EVLAB toolbox version 17 or above and try again'); end
                end
                if strcmp(pwd,fileparts(which(mfilename))),
                    fprintf('warning: adding current directory %s to Matlab path; please include the FL software folder in your Matlab path and cd to a different folder to avoid this warning in the future\n',pwd);
                    addpath(pwd);
                end
                done_init=path;
            end
            evlab17(STEPS,varargin{:});
        case {'versioncheckreset'}
            fprintf('FL will now continue performing software-version number/release/folder checks when run\n');
            done_init=[];
        case {'versioncheckskip'}
            fprintf('warning: FL will now stop performing software-version number/release/folder checks when run\n');
            fprintf('Please use "fl_process versioncheckreset", or start a new Matlab session, to continue performing software-version checks\n');
            done_init=0;
        case 'utils_rootfolder'
            if nargin>1&&~isempty(varargin{1}), rootfolder=varargin{1}; end
            if isempty(rootfolder), varargout={fullfile(fileparts(which(mfilename)),'REPOSITORY')};
            else varargout={rootfolder};
            end
        case 'utils_setpermissions'
            fl_process init silent;
            evlab17_setpermissions;
            evlab17_fixpermissions(varargin{:});            
        case {'qa.plots','qaplots'}
            fl_process init silent;
            conn_qaplotsexplore(varargin{:},'flfolders');
        case 'open'
            fl_process init silent;
            evlab17_module('load',varargin{1});
            conn;
            conn('load',conn_prepend('',varargin{1},'.mat'));
            conn gui_setup;
        case {'model.plots','modelplots'}
            fl_process init silent;
            evlab17_modelplots(varargin{:});
        case {'parallel.report','report'}
            fl_process init silent;
            evlab17_module('load',varargin{1});
            conn_jobmanager report;
        case {'parallel.report.gui','report.gui'}
            fl_process init silent;
            evlab17_module('load',varargin{1});
            conn_jobmanager;
        case {'parallel.cancel','cancel'}
            fl_process init silent;
            evlab17_module('load',varargin{1});
            conn_jobmanager cancel;
        case {'parallel.delete','delete'}
            fl_process init silent;
            evlab17_module('load',varargin{1});
            conn_jobmanager delete;
        otherwise
            if ~isempty(which(sprintf('evlab17_%s',STEPS))),
                fh=eval(sprintf('@evlab17_%s',STEPS));
                if ~nargout, feval(fh,varargin{:});
                else [varargout{1:nargout}]=feval(fh,varargin{:});
                end
            elseif ~isempty(which(sprintf('conn_%s',STEPS))),
                fh=eval(sprintf('@conn_%s',STEPS));
                if ~nargout, feval(fh,varargin{:});
                else [varargout{1:nargout}]=feval(fh,varargin{:});
                end
            else
                disp(sprintf('unrecognized option %s',STEPS));
            end
    end
    return
end

fl_process init silent;

STEPS=cellstr(STEPS);
if rem(numel(varargin),2), OPTS=reshape([{'subject_id'},varargin],2,[]);
else OPTS=reshape(varargin,2,[]);
end

OVERWRITE=true;
SUBJECT_ID={};
SUBJECT_INFO={};
PIPELINE_INFO={};
DESIGN_INFO={};
OUTPUT_FOLDER=fl_process('utils_rootfolder');
OUTPUT_FOLDER_subject='';
OPTS_remove=true(1,size(OPTS,2));
PARALLEL=false; %any(ismember(lower(STEPS),{'all','preproc','preproc.append'}));
LOCALCOPY=true;
QAPLOTS=false;
QAPLIST=[];
QASET=[];
IMMEDIATERETURN=numel(STEPS)==1&&all(ismember(lower(STEPS),{'preproc','preproc.append','qa.create','qacreate'}));
DONOTEXPAND=false;
for nopts=1:size(OPTS,2), 
    switch(lower(OPTS{1,nopts}))
        case 'subject_info',        SUBJECT_INFO=OPTS{2,nopts}; 
        case 'subject_id',          SUBJECT_ID=OPTS{2,nopts};
        case 'pipeline_info',       PIPELINE_INFO=OPTS{2,nopts}; 
        case 'design_info',         DESIGN_INFO=OPTS{2,nopts}; 
        case 'root_folder',         OUTPUT_FOLDER=OPTS{2,nopts};
        case 'overwrite',           OVERWRITE=OPTS{2,nopts};
        case 'parallel',            PARALLEL=OPTS{2,nopts};
        case 'qa_plist',            QAPLIST=OPTS{2,nopts};
        case 'qa_set',              QASET=OPTS{2,nopts};
        case 'qa_plots',            QAPLOTS=OPTS{2,nopts};
        case 'localcopy',           LOCALCOPY=OPTS{2,nopts};
        case 'immediatereturn',     IMMEDIATERETURN=OPTS{2,nopts};
        case 'donotexpand',         DONOTEXPAND=OPTS{2,nopts};   % if single subject_id string and multipel subject_info files: donotexpand=0: expands SUBJECT_ID to SUBJECT_ID/1 SUBJECT_ID/2 ...; donotexpand=1: expands SUBJECT_ID to SUBJECT_ID,1 SUBJECT_ID,2 ...
        otherwise,                  OPTS_remove(nopts)=false;
    end
end
OPTS(:,find(OPTS_remove))=[];
OPTS=reshape(OPTS,1,[]);
if ischar(PARALLEL), PARALLEL=str2num(PARALLEL); end
if ischar(OVERWRITE), OVERWRITE=str2num(OVERWRITE); end
if ischar(LOCALCOPY), LOCALCOPY=str2num(LOCALCOPY); end
if ischar(QAPLOTS), QAPLOTS=str2num(QAPLOTS); end
if ischar(QAPLIST)&&~isempty(str2num(QAPLIST)), QAPLIST=str2num(QAPLIST); end
if ischar(IMMEDIATERETURN), IMMEDIATERETURN=str2num(IMMEDIATERETURN); end
if ischar(DONOTEXPAND), DONOTEXPAND=str2num(DONOTEXPAND); end
if PARALLEL, OPTS=[OPTS,{'parallel.N',1,'parallel.immediatereturn',IMMEDIATERETURN}]; end
if LOCALCOPY, OPTS=[OPTS,{'localcopy',1,'localcopy_reduce',1}]; end
OPTS=[OPTS,{'qa_plots',QAPLOTS}];
if ~isempty(QAPLIST), OPTS=[OPTS,{'qa_plist',QAPLIST}]; end
if ~isempty(QASET), OPTS=[OPTS,{'qa_set',QASET}]; end

if ischar(SUBJECT_ID)&&any(ismember(SUBJECT_ID,'*?'))
    temp=conn_dir(SUBJECT_ID,'-R','-cell','-dir');
    assert(~isempty(temp),'No match to %s',SUBJECT_ID);
    if isempty(temp)
        fprintf('no ID files match found (%s). Attempting to derive ID from subject data information filenames\n',SUBJECT_ID);
    else
        SUBJECT_ID=conn_prepend('',temp,'.mat');
        fprintf('ID files match: %s\n',sprintf('%s ',SUBJECT_ID{:}));
    end
end
if ischar(SUBJECT_INFO)&&any(ismember(SUBJECT_INFO,'*?'))
    temp=conn_dir(SUBJECT_INFO,'-R','-cell');
    assert(~isempty(temp),'No match to %s',SUBJECT_INFO);
    SUBJECT_INFO=temp;
    fprintf('CFG files match: %s\n',sprintf('%s ',SUBJECT_INFO{:}));
end
if isempty(SUBJECT_ID)||(ischar(SUBJECT_ID)&&any(ismember(SUBJECT_ID,'%*')))
    assert(~isempty(SUBJECT_INFO),'No subject data information files selected');
    [nill,temp,nill]=cellfun(@fileparts,cellstr(SUBJECT_INFO),'uni',0);
    if ~isempty(regexp(SUBJECT_ID,'\%s')), SUBJECT_ID=cellfun(@(x)sprintf(SUBJECT_ID,x),temp,'uni',0); % e.g. SUBJECT_ID='sub-%s';
    elseif ~isempty(regexp(SUBJECT_ID,'\*')), SUBJECT_ID=cellfun(@(x)regexprep(SUBJECT_ID,'[^\/\\]*\*[^\/\\]*',x),temp,'uni',0); % e.g. SUBJECT_ID='sub-*';
    elseif ~isempty(regexp(SUBJECT_ID,'\%\d+d')), SUBJECT_ID=cellfun(@(x)sprintf(SUBJECT_ID,str2double(regexp(x,'\d+','match','once'))),temp,'uni',0); % e.g. SUBJECT_ID='sub-%04d';
    else SUBJECT_ID=temp;
    end
    fprintf('Subject IDs automatically set to %s\n',sprintf('%s ',SUBJECT_ID{:}));
end
assert(~isempty(SUBJECT_ID),'No subject ID entered');
if ~iscell(SUBJECT_ID), SUBJECT_ID=cellstr(SUBJECT_ID); end
if ~iscell(SUBJECT_INFO),   if ischar(SUBJECT_INFO), SUBJECT_INFO=cellstr(SUBJECT_INFO); else SUBJECT_INFO={SUBJECT_INFO}; end; end
if ~iscell(PIPELINE_INFO),  if ischar(PIPELINE_INFO), PIPELINE_INFO=cellstr(PIPELINE_INFO); else PIPELINE_INFO={PIPELINE_INFO}; end; end
if ~iscell(DESIGN_INFO),    if ischar(DESIGN_INFO), DESIGN_INFO=cellstr(DESIGN_INFO); else DESIGN_INFO={DESIGN_INFO}; end; end
assert(~any(ismember(lower(STEPS),'preproc.append'))||isempty(SUBJECT_INFO),'Incorrect usage: PREPROC.APPEND cannot be combined with subject_info field; use PREPROC to import&preprocess a new dataset, or remove subject_info field to preprocess an already imported dataset');
assert(~any(ismember(lower(STEPS),'preproc.import'))||isempty(PIPELINE_INFO),'Incorrect usage: PREPROC.IMPORT cannot be combined with pipeline_info field; use PREPROC to import&preprocess a new dataset, or remove pipeline_info field to import a new dataset');
if any(ismember(lower(STEPS),{'preproc','preproc.import','all'}))
    if numel(SUBJECT_ID)==1&&numel(SUBJECT_INFO)>1, 
        if DONOTEXPAND, SUBJECT_ID=arrayfun(@(x)sprintf('%s,%d',SUBJECT_ID{1},x),1:numel(SUBJECT_INFO),'uni',0); % keep in single project
            fprintf('Warning: %d ID files with %d CFG files. Found ''donotexpand''=%d setting. Creating %d-subjects single project %s\n',numel(SUBJECT_ID),numel(SUBJECT_INFO),DONOTEXPAND,numel(SUBJECT_INFO));
        else
            [nill,subject_info_names,nill]=cellfun(@fileparts,SUBJECT_INFO,'uni',0);
            SUBJECT_ID=cellfun(@(x)sprintf('%s/%s',conn_prepend('',SUBJECT_ID{1},''),x),subject_info_names,'uni',0);                 % expand to multiple projects
            fprintf('Warning: %d ID files with %d CFG files. Found ''donotexpand''=%d setting. Expanding to multiple projects %s\n',numel(SUBJECT_ID),numel(SUBJECT_INFO),DONOTEXPAND,sprintf('%s ',SUBJECT_ID{:}));
        end
    end
    assert(numel(SUBJECT_INFO)==numel(SUBJECT_ID),'Mismatched number of subjects (%d subject information files, %d subject ids)',numel(SUBJECT_INFO),numel(SUBJECT_ID));
end
if any(ismember(lower(STEPS),{'preproc','preproc.append','all'}))
    assert(numel(SUBJECT_ID)==numel(PIPELINE_INFO)|numel(PIPELINE_INFO)<=1,'Mismatched number of preprocessing information files (%d subject ids, %d preprocessing information files)',numel(SUBJECT_ID),numel(PIPELINE_INFO));
    if numel(PIPELINE_INFO)==1, PIPELINE_INFO=repmat(PIPELINE_INFO,1,numel(SUBJECT_ID)); end
end
if any(ismember(lower(STEPS),{'model','all'}))
    assert(numel(SUBJECT_ID)==numel(DESIGN_INFO)|numel(DESIGN_INFO)<=1,'Mismatched number of design information files (%d subject ids, %d design information files)',numel(SUBJECT_ID),numel(DESIGN_INFO));
    if numel(DESIGN_INFO)==1, DESIGN_INFO=repmat(DESIGN_INFO,1,numel(SUBJECT_ID)); end
end
if isempty(OUTPUT_FOLDER), OUTPUT_FOLDER=fl_process('utils_rootfolder');; end
[OUTPUT_FOLDER_subject,SUBJECT_ID]=cellfun(@fileparts,SUBJECT_ID,'uni',0);

pwd0=pwd;
filesout={};


%if ~isempty(dir('logfile.mat')), load('logfile.mat','filesout'); end
previous_subject_id={};
for nsub=1:numel(SUBJECT_ID)
    if numel(OUTPUT_FOLDER_subject)>=nsub&&~isempty(OUTPUT_FOLDER_subject{nsub}), 
        if ~isempty(regexp(OUTPUT_FOLDER_subject{nsub},'^(\.)')), pwd1=fullfile(pwd0,OUTPUT_FOLDER_subject{nsub}); % subject_id contains relative path
        elseif ~isempty(regexp(OUTPUT_FOLDER_subject{nsub},'^(\/)')), pwd1=OUTPUT_FOLDER_subject{nsub}; % subject_id contains full path
        else pwd1=fullfile(OUTPUT_FOLDER,OUTPUT_FOLDER_subject{nsub});             % subject_id contains partial path
        end
    else pwd1=OUTPUT_FOLDER;                                                       % subject_id contains no path info
    end
    [nill,nill]=mkdir(pwd1);
    subject_id=conn_prepend('',SUBJECT_ID{nsub},'');
    assert(~ismember(subject_id,previous_subject_id),'Found repeated subject_id %s. Analysis stopped',subject_id);
    previous_subject_id=[previous_subject_id {subject_id}];
    thisOPTS=OPTS;
    if ~isempty(regexp(subject_id,',\d+$'))
        dataset_subject=str2num(char(regexp(subject_id,',(\d+)$','tokens','once')));
        thisOPTS=[thisOPTS, {'subjects',dataset_subject}];
        subject_id=regexprep(subject_id,',\d+$','');
    end
    dataset=conn_prepend('',fullfile(pwd1,subject_id),'.mat');
    for ntin=1:2
        if ntin==1,  tin=PIPELINE_INFO;
        else         tin=SUBJECT_INFO;
        end
        if numel(tin)>=nsub&&~isempty(tin{nsub})&&ischar(tin{nsub})
            if ~isempty(regexp(tin{nsub},'^(\.)')), tout=fullfile(pwd0,tin{nsub});
            elseif isempty(regexp(tin{nsub},'^(\/)')), tout=fullfile(OUTPUT_FOLDER,tin{nsub});
            else tout=tin{nsub};
            end
            if ntin==1,  PIPELINE_INFO{nsub}=tout;
            else         SUBJECT_INFO{nsub}=tout;
            end
        end
    end
    
    if any(ismember(lower(STEPS),{'all','preproc'}))
        %try
            if OVERWRITE||isempty(dir(dataset))
                if isempty(PIPELINE_INFO), pipeline={};
                else
                    pipeline=PIPELINE_INFO{nsub};
                    if ~isempty(pipeline)&&ischar(pipeline)&&isempty(fileparts(pipeline))&&~isempty(dir(fullfile(fileparts(which(mfilename)),pipeline))), pipeline=fullfile(fileparts(which(mfilename)),pipeline); end
                    pipeline={pipeline};
                end
                filesout{nsub}=dataset;
                if ischar(SUBJECT_INFO{nsub})&&~isempty(regexp(SUBJECT_INFO{nsub},'\.mat$')) % if subject info is another dataset
                    dataset0=SUBJECT_INFO{nsub};
                    evlab17_module('load',dataset0);
                    pwd2=fileparts(dataset);
                    if ~isdir(pwd2), [nill,nill]=mkdir(pwd2); end
                    evlab17_module('save',dataset); 
                    evlab17_module('update'); % import
                    if LOCALCOPY, 
                        conn_importvol2bids(true); 
                        conn save;
                    end
                    evlab17_run_preproc(... % & append
                        pipeline{:},...
                        [],...
                        thisOPTS{:},...
                        'dataset',dataset);
                else % if subject info is a cfg file or structure
                    evlab17_run_preproc(...
                        SUBJECT_INFO{nsub},...
                        pipeline{:},...
                        [],...
                        thisOPTS{:},...
                        'dataset',dataset);
                end
            else fprintf('Skipping preprocessing subject %s\n',SUBJECT_ID{nsub}); cd(pwd0);
            end
        %catch me
        %    fprintf('----------------------------\n\nERROR in %s\nSKIPPING SUBJECT\n%s\n\n',dataset,me.message); cd(pwd0)
        %end
    end
       
    if any(ismember(lower(STEPS),{'preproc.import'}))
        %try
            if OVERWRITE||isempty(dir(dataset))
                filesout{nsub}=evlab17_run_preproc(...
                    SUBJECT_INFO{nsub},...
                    [],...
                    thisOPTS{:},...
                    'dataset',dataset);
            else fprintf('Skipping import subject %s\n',SUBJECT_ID{nsub}); cd(pwd0);
            end
        %catch me
        %    fprintf('----------------------------\n\nERROR in %s\nSKIPPING SUBJECT\n%s\n\n',dataset,me.message); cd(pwd0)
        %end
    end
    
    
    if any(ismember(lower(STEPS),{'preproc.append'}))
        %try
            if OVERWRITE||isempty(dir(dataset))
                if isempty(PIPELINE_INFO), pipeline={};
                else
                    pipeline=PIPELINE_INFO{nsub};
                    if ~isempty(pipeline)&&ischar(pipeline)&&isempty(fileparts(pipeline))&&~isempty(dir(fullfile(fileparts(which(mfilename)),pipeline))), pipeline=fullfile(fileparts(which(mfilename)),pipeline); end
                    pipeline={pipeline};
                end
                filesout{nsub}=dataset;
                evlab17_run_preproc(...
                    pipeline{:},...
                    [],...
                    thisOPTS{:},...
                    'dataset',dataset);
            else fprintf('Skipping preprocessing subject %s\n',SUBJECT_ID{nsub}); cd(pwd0);
            end
        %catch me
        %    fprintf('----------------------------\n\nERROR in %s\nSKIPPING SUBJECT\n%s\n\n',dataset,me.message); cd(pwd0)
        %end
    end
       
    if any(ismember(lower(STEPS),{'all','model'}))
        %try
            if OVERWRITE||isempty(dir(dataset))
                if isempty(DESIGN_INFO), pipeline={};
                else
                    pipeline=DESIGN_INFO{nsub};
                    if ~isempty(pipeline)&&ischar(pipeline)&&isempty(fileparts(pipeline))&&~isempty(dir(fullfile(fileparts(which(mfilename)),pipeline))), pipeline=fullfile(fileparts(which(mfilename)),pipeline); end
                    pipeline={pipeline};
                end
                filesout{nsub}=dataset;
                if any(ismember(lower(STEPS),{'all'}))||isempty(SUBJECT_INFO), 
                    evlab17_module('load',dataset);
                    options=evlab17_module('getinfo','design');
                    assert(isfield(options,'design')||isfield(options,'files'),'design information not found in original dataset %s',dataset);
                    evlab17_run_model(...
                        struct('model_session',0),... % note: changes default from evlab17_run_model
                        options,...
                        pipeline{:},...
                        [],...
                        thisOPTS{:},...
                        'dataset',dataset);
                else %note: this allows design info to be temporarily modified when running model step
                    fprintf('warning: disregarding original dataset design info, using information in %s instead\n',SUBJECT_INFO{nsub});
                    evlab17_run_model(...
                        struct('model_session',0),... % note: changes default from evlab17_run_model
                        SUBJECT_INFO{nsub},...
                        pipeline{:},...
                        [],...
                        thisOPTS{:},...
                        'dataset',dataset);
                end
            else fprintf('Skipping first-level analysis subject %s\n',SUBJECT_ID{nsub}); cd(pwd0);
            end
        %catch me
        %    fprintf('----------------------------\n\nERROR in %s\nSKIPPING SUBJECT\n%s\n\n',dataset,me.message); cd(pwd0)
        %end
    end
    
    if any(ismember(lower(STEPS),{'all','qa.create','qacreate'}))
        if OVERWRITE||isempty(dir(dataset))
            idx=2*find(ismember(thisOPTS(1:2:end-1),{'subjects','qa_plist','qa_set','parallel.N','parallel.immediatereturn','parallel.profile','qa_parallel','qa_profile'}))-1;
            opts=thisOPTS([idx(:) idx(:)+1]');
            filesout{nsub}=evlab17_run_qa(...
                [],...
                'dataset',dataset,...
                opts{:});
                %'qa_parallel',0,...
        else fprintf('Skipping Quality Control plots subject %s\n',SUBJECT_ID{nsub});
        end
    end
    
    if any(ismember(lower(STEPS),{'qa.plot','qaplot','qa'}))
        if conn_existfile(dataset)
            evlab17_module('load',dataset);
            conn_qaplotsexplore;
        elseif isdir(fullfile(pwd1,subject_id))
            conn_qaplotsexplore(fullfile(pwd1,subject_id))
        end
    end
end
varargout{1}=filesout;
%cd(pwd0);

end



