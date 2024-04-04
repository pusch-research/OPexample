ccc
warning('off','MATLAB:datetime:FormatConflict_mM')

% define paths
paths.specs=fullfile('OpenFAST-DTU10MW-model','specs','specsDTU10MW.m');
paths.CpCqCt=fullfile('OpenFAST-DTU10MW-model','specs','Cp_Ct_Cq.DTU10MW_Nautilus.txt');
paths.FSTtmplInputFile=fullfile('OpenFAST-DTU10MW-model','template','Model.fst');
paths.FSTnewInputFileDir='C:\tmp\OpenFAST-DTU10MW';
paths.FSTbin='OpenFAST-bin-v3-WINx64';
paths.results='results';
paths.resultsTmpDir=paths.FSTnewInputFileDir;
paths.parforProgressFile=fullfile(getPath(paths.FSTnewInputFileDir),'parfor_progress.txt');

% run init
initSim


% open simulink simulation model
% open(simOptTmpl.SLmodelName);


%% overwrite simOpt default parameters

% simOptTmpl.settlingAbsThresh=0.02;
% simOptTmpl.FSTdataTmpl.DT=0.025; % needs to be 0.025


%% run simulations (to generate raw data)

% define meta data
metaData=struct();
metaData.dateStr=char(datetime('now','Format','yy-mm-dd_HHmm'));
metaData.modelName='DTU10MW';
metaData.info='default';


% outer loop for windSpeed and plattform pitch
windSpeed_mDs_arr=10; %unique([4 5 6 7 8 8.5 9 9.5 10 10.5 11:0.1:12 12:0.5:15 15:24]);
ptfmPitch_deg_arr=unique([0 2]); % [deg]
n_windSpeed=numel(windSpeed_mDs_arr);
n_ptfmPitch=numel(ptfmPitch_deg_arr);
n_outer=n_windSpeed*n_ptfmPitch;
% inner loop for GenTq and BldPitchC
n_GenTqSmpl=3;%
n_BldPitchCSmpl=3; %


% start loop
tic
parfor_progress(n_outer,getPath(paths.parforProgressFile));
res_arr=struct('outAvg',[],'settlingTime',[],'isValid',[],'windSpeed',[],'GenTq',[],'BldPitchC',[],'PtfmPitch',[]);
res_arr=repmat(res_arr,n_windSpeed,n_ptfmPitch);

% Note: here you can select if you want to run a FOR-loop or a PARFOR-loop.
% To do so, un/comment the according lines!

% % ##### PARFOR loop
% myParPool=gcp('nocreate'); 
% if isempty(myParPool)
%     myParPool=parpool(min(n_outer,floor(str2double(getenv('NUMBER_OF_PROCESSORS'))*0.9)));
% end
% parfor i_outer=1:n_outer
  

% ##### FOR loop
for i_outer=1:n_outer
    

    % re-initialize (recommended for parfor)
    simOpt=simOptTmpl; 

    % command line output
    [i_wind,i_ptfmPitch]=ind2sub([n_windSpeed n_ptfmPitch],i_outer);
    windSpeed_mDs_act=windSpeed_mDs_arr(i_wind);
    ptfmPitch_deg_act=ptfmPitch_deg_arr(i_ptfmPitch);
    disp(['####################### START i=' num2str(i_outer) ...
          ' (windSpeed=' num2str(windSpeed_mDs_act) 'm/s, ptfmPitch=' num2str(ptfmPitch_deg_act) 'deg) '...
          '#########################  (' char(datetime('now','Format','HH:mm')) ')']) 
    
    % compute upper and lower bounds for gridding for given wind speed
    % based on reference (= theoretical optimal) values 
    % this involves quite some trial and error and is tailored to the DTU10MW wind turbine
    RotSpeedRef_rpm=interp1(specs.schedule.v_mDs_arr,specs.schedule.wr_radDs_arr,windSpeed_mDs_act)*radDs2rpm; 
    GenTqRef_Nm=interp1(specs.schedule.v_mDs_arr,specs.schedule.GenTq_Nm_arr,windSpeed_mDs_act); 
    BldPitchCRef_rad=interp1(specs.schedule.v_mDs_arr,specs.schedule.BldPitch_deg_arr,windSpeed_mDs_act)*deg2rad(1);
    if windSpeed_mDs_act<specs.RtdWndSpeed
        % below rated
        GenTq_Nm_lim=[0.9 1.2]*GenTqRef_Nm; % 0.9 1.1 [0.95 1.3] [0.95 1.2]
        BldPitchC_rad_lim=[-2 1.2]*BldPitchCRef_rad; % [-3 1.2]
    elseif windSpeed_mDs_act<specs.RtdWndSpeed+1
        % region 2.5
        GenTq_Nm_lim=[0.8 1.2]*GenTqRef_Nm;
        BldPitchC_rad_lim=deg2rad([-0.2+specs.schedule.BldPitch_deg_arr(1) ...
                                   max(rad2deg(BldPitchCRef_rad)+0.2,specs.schedule.BldPitch_deg_arr(1)+3)]);
    else
        % above rated
        GenTq_Nm_lim=[0.8 1.2]*GenTqRef_Nm;
        BldPitchC_rad_lim=[0.9*interp1(specs.schedule.v_mDs_arr,specs.schedule.BldPitch_deg_arr,windSpeed_mDs_act-0.5,'linear','extrap')*deg2rad(1) ...
                           1.2*interp1(specs.schedule.v_mDs_arr,specs.schedule.BldPitch_deg_arr,windSpeed_mDs_act+1.5,'linear','extrap')*deg2rad(1)];
    end

    % define simulink input (grid points for GenTq and BldPitchC around theoretical optimal values)
    SLin=struct();
    if n_GenTqSmpl>1
        n_tmp=min(n_GenTqSmpl-1,max(2,round((n_GenTqSmpl+1)*(GenTqRef_Nm-GenTq_Nm_lim(1))/range(GenTq_Nm_lim))));
        GenTq_Nm_arr=unique([linspace(GenTq_Nm_lim(1),GenTqRef_Nm,n_tmp)...
                             linspace(GenTqRef_Nm,GenTq_Nm_lim(2),(n_GenTqSmpl+1)-n_tmp)]);
    else
        GenTq_Nm_arr=GenTqRef_Nm;
    end
    if n_BldPitchCSmpl>1
        n_tmp=min(n_BldPitchCSmpl-1,max(2,round((n_BldPitchCSmpl+1)*(BldPitchCRef_rad-BldPitchC_rad_lim(1))/range(BldPitchC_rad_lim))));
        BldPitchC_rad_arr=unique([linspace(BldPitchC_rad_lim(1),BldPitchCRef_rad,n_tmp)...
                                  linspace(BldPitchCRef_rad,BldPitchC_rad_lim(2),(n_BldPitchCSmpl+1)-n_tmp)]);
    else
        BldPitchC_rad_arr=BldPitchCRef_rad;
    end
    [SLin.BldPitchC,SLin.GenTq]=ndgrid(BldPitchC_rad_arr,GenTq_Nm_arr); % note; first vary BldPitch, then increase GenTq
    SLin.n_step=numel(SLin.BldPitchC); % number of parameter combinations; assume that x_scaled is a row vector (or array of row vectors)
    SLin.extra=zeros(43,1); % extra inputs (e.g. cable tensions, flaps, etc.)


    % update simOpt (and generate new OpenFAST input files)
    FSTdata=simOpt.FSTdataTmpl;
    FSTdata.InflowFile.HWindSpeed=windSpeed_mDs_act;
    FSTdata.ElastoFile.RotSpeed=RotSpeedRef_rpm;  % overwrite initial value in ElastoDyn
    FSTdata.ElastoFile.PtfmPitch=ptfmPitch_deg_act;
    FSTdata.TMax=simOpt.SLmaxTime;
    simOpt.FSTnewModelSuffix=num2str(['_W' num2str(windSpeed_mDs_act) 'PP' num2str(ptfmPitch_deg_act)]);
    [simOpt.FSTnewInputFile,tmpFSTfileName_arr]=generateFASTinFiles(...
        simOpt.FSTtmplInputFile,...
        simOpt.FSTnewInputFileDir,...
        simOpt.FSTnewModelSuffix,FSTdata);  
    simOpt.settlingAbsThresh=RotSpeedRef_rpm*0.0007;

    % check if resTmpFile already exists and load it in case (useful, e.g., when simulation has crashed)
    resTmpFileName=fullfile(getPath(paths.resultsTmpDir),['TMP' num2str(i_outer) simOpt.FSTnewModelSuffix '.mat']);
    if exist(resTmpFileName,'file')
        
        % resTmpFile exists and is loaded
        try
            res_act=loadvar(resTmpFileName,'res_act');
            disp(['> ' resTmpFileName ' loaded...'])
        catch
            delete(resTmpFileName)
            warning(['> ERROR loading ' resTmpFileName '. File deleted.'])
        end

    end
    if ~exist(resTmpFileName,'file')
        
        % resTmpFile does not exist - run simulation
        simout=doSim(simOpt,SLin);
    
        % store results
        res_act=struct();
        res_act.outAvg=simout.outAvg;
        res_act.settlingTime=simout.settlingTime;
        res_act.isValid=simout.isValid;
        res_act.windSpeed=windSpeed_mDs_act; % m/s
        res_act.GenTq=SLin.GenTq; % Nm
        res_act.BldPitchC=SLin.BldPitchC; % rad
        res_act.PtfmPitch=ptfmPitch_deg_act; % deg
    
        % tmp save results
        parsave(resTmpFileName,'res_act',res_act,'simOpt',simOpt)
        disp(['> ' resTmpFileName ' saved..'])

    end
    
    % store in array
    res_arr(i_outer)=res_act;

    % clean-up and show progress
    delete(tmpFSTfileName_arr{:}); % delete temporary OpenFAST input files
    percentDone=parfor_progress([],getPath(paths.parforProgressFile));
    disp(['####################### FINISH i=' num2str(i_outer) ...
          ' (windSpeed=' num2str(windSpeed_mDs_act) 'm/s,ptfmPitch=' num2str(ptfmPitch_deg_act) 'deg) '...
          '######################### (' char(datetime('now','Format','HH:mm')) ') ' num2str(percentDone) '%']) 

end
elapsedTime=toc;
parfor_progress(0,getPath(paths.parforProgressFile)); % reset progress bar (delete .txt file)

% % save res_arr for later usage (if desired)
% resTmpFileName=fullfile(getPath(paths.resultsTmpDir),'TMPres_arr.mat');
% save(resTmpFileName,'res_arr','simOptTmpl','specs'); disp(['> ' resTmpFileName ' saved...']);
% load(resTmpFileName)


%% post-process simulation data (for efficient data analysis)

% data_arr is a [n_wind 1] cell array with each cell being a 
% struct with fields
%   isValid     numeric array [n_BldPitch n_GenTq n_ptfmPitch]
%   [OUTNAME]   numeric array [n_BldPitch n_GenTq n_ptfmPitch]
%   windSpeed   scalar in m/s
%   dim         struct with information about the grid dimensions

[n_wind,n_ptfmPitch]=size(res_arr);
n_out=numel(simOptTmpl.FSToutNameArr);
data_arr=cell(n_wind,1);
for i_wind=1:n_wind

    [n_BldPitch,n_GenTq]=size(res_arr(i_wind,1).BldPitchC);
    data_arr{i_wind}.isValid=nan(n_BldPitch,n_GenTq,n_ptfmPitch);
    for i_ptfmPitch=1:n_ptfmPitch
        for i_out=1:n_out
            outName_act=simOptTmpl.FSToutNameArr{i_out};
            data_arr{i_wind}.(outName_act)(:,:,i_ptfmPitch)=reshape(res_arr(i_wind,i_ptfmPitch).outAvg(:,i_out),n_BldPitch,n_GenTq);
        end
        data_arr{i_wind}.isValid(:,:,i_ptfmPitch)=reshape(res_arr(i_wind,i_ptfmPitch).isValid,n_BldPitch,n_GenTq);
    end
    data_arr{i_wind}.windSpeed=res_arr(i_wind).windSpeed;
    data_arr{i_wind}.grid.dimNames={'BldPitchC','GenTq','PtfmPitch'};
    data_arr{i_wind}.grid.dimSmplPts={ res_arr(i_wind,1).BldPitchC(:,1)'*180/pi ...
                                     res_arr(i_wind,1).GenTq(1,:) ...
                                     [res_arr(i_wind,:).PtfmPitch]};
    data_arr{i_wind}.grid.BldPitchC=res_arr(i_wind,1).BldPitchC(:,1)'*180/pi;
    data_arr{i_wind}.grid.GenTq=res_arr(i_wind,1).GenTq(1,:);
    data_arr{i_wind}.grid.PtfmPitch=[res_arr(i_wind,:).PtfmPitch];
    data_arr{i_wind}.grid.dimUnits={'deg' 'Nm' 'deg'};
    for ii=1:numel(data_arr{i_wind}.grid.dimNames)
        dimName_act=data_arr{i_wind}.grid.dimNames{ii};
        data_arr{i_wind}.grid.dimIdx.(dimName_act)=ii;
    end

    % HACK: compute TSR manually
    data_arr{i_wind}.TSR=data_arr{i_wind}.RotSpeed*rpm2radDs*specs.R/data_arr{i_wind}.windSpeed;


end
n_grid=numel(data_arr{1}.grid.dimNames);



%% save and clean-up

% save results
if true
    metaData.elapsedTime=elapsedTime;
    resFileName=fullfile(getPath(paths.results),[metaData.dateStr '_' metaData.info '.mat']);
    save(resFileName,'data_arr','simOptTmpl','paths','specs','metaData')
    disp(['> ' resFileName ' saved..'])
end

% empty tmp result files folder
if true
    tmpDirInfo = dir(getPath(paths.resultsTmpDir));
    tmpDirInfo([tmpDirInfo.isdir]) = [];   % skip directories
    filenames = fullfile(getPath(paths.resultsTmpDir), {tmpDirInfo.name});
    delete( filenames{:} )
    disp('> tmp result files deleted..')
end

disp(['done in ' seconds2human(elapsedTime) ' (' char(datetime('now','Format','HH:mm')) ')']);


%% do data analysis (missing)

