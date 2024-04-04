ccc
load('results\22-08-06_1625_fixedBottom.mat')

% re-assign
windSpeed_arr=cellfun(@(x) x.windSpeed,data_arr);
n_wind=numel(windSpeed_arr);

%% modify data
% HACK: replace invalid values with nearest valid values
% if you don't do this you run into various problems during optimization
outName_arr={'GenPwr' 'GenSpeed' 'RotThrust' 'TSR' 'RotSpeed' 'TwrBsMyt' }; % simOptTmpl.FSToutNameArr;
for i_wind=1:n_wind
    for i_out=1:numel(outName_arr)
        outName_act=outName_arr{i_out};
        out_act=data_arr{i_wind}.(outName_act);
        out_act(~data_arr{i_wind}.isValid)=nan;
        data_arr{i_wind}.(outName_act)=fillmissing(out_act,'nearest',2,'EndValues','nearest');
    end
    data_arr{i_wind}.isValid=true(size(data_arr{i_wind}.isValid)); % HACK
end

% add objective function (to be maximized) for making tower base loads zero
for i_wind=1:n_wind
    data_arr{i_wind}.TwrBsMytZERO=-abs(data_arr{i_wind}.TwrBsMyt);
end


%% compute optimal schedule

wind_iArr=1:n_wind; % only below rated

% options
interpMethod='makima'; % spline is better than makima
fminconOpt=[];
% select parameters
pVal_arr={data_arr{1}.grid.PtfmPitch}; % to be optimized as well
% select output channels for schedule
schedOutName_arr={'GenSpeed' 'RotThrust' 'TSR' 'RotSpeed' 'TwrBsMyt' 'RtTSR'};

% optimize all: BldPitchC, GenTq, and PtfmPitch
disp('> compute optimal schedule (all)..')
ctrlOpt=struct();
ctrlOpt.belowRtd=struct('GenTq','opt','BldPitchC','opt');
ctrlOpt.belowRtd.cnstr=[];
ctrlOpt.aboveRtd=struct('GenTq','rtd','BldPitchC','optTwrBsMytZERO'); % optTwrBsMytZERO,opt
ctrlOpt.aboveRtd.cnstr=struct('OutName',{{'PtfmPitch'}},'ub',inf,'lb',nan);
[scheduleOpt,rtdOpt]=scheduleInterpGridData(data_arr(wind_iArr),schedOutName_arr,pVal_arr,specs,ctrlOpt,fminconOpt,interpMethod);

% % optimize BldPitchC and PtfmPitch + Kw2 law (DOES NOT WORK)
% disp('> compute optimal schedule (BldPitchC only) ..')
% ctrlBP=ctrlOpt;
% ctrlBP.belowRtd=struct('GenTq','Kw2','BldPitchC','opt');
% [scheduleBP,rtdBP]=scheduleInterpGridData(data_arr(wind_iArr),schedOutName_arr,pVal_arr,specs,ctrlBP,fminconOpt,interpMethod);

% optimize GenTq and PtfmPitch + fine pitch
disp('> compute optimal schedule (GentTq only) ..')
ctrlGT=ctrlOpt;
ctrlGT.belowRtd=struct('GenTq','opt','BldPitchC','finePitch');
[scheduleGT,rtdGT]=scheduleInterpGridData(data_arr(wind_iArr),schedOutName_arr,pVal_arr,specs,ctrlGT,fminconOpt,interpMethod);

% optimize PtfmPitch + Kw2 law + fine pitch
disp('> compute optimal schedule (PtfmPitch only) ..')
ctrlPP=ctrlOpt;
ctrlPP.belowRtd=struct('GenTq','Kw2','BldPitchC','finePitch');
[schedulePP,rtdPP]=scheduleInterpGridData(data_arr(wind_iArr),schedOutName_arr,pVal_arr,specs,ctrlPP,fminconOpt,interpMethod);


%% plot BldPitch-GenTq CONTOURF for single wind speed and platform pitch

% close all

xName='BldPitchC';
yName='GenTq';
zName='GenPwr'; % GenPwr, TwrBsMyt

PtfmPitch_deg=0; % may be overwritten inside loop
infoStr=['PtfmPitch=' num2str(PtfmPitch_deg) 'deg'];



for i_wind=10 %1:n_wind

    data_act=data_arr{i_wind};

    % overwrite ptfmPitch by value from david
%     PtfmPitch_deg=interp1(specs.refSchedule.windSpeed,specs.refSchedule.PtfmPitch,data_act.windSpeed); % ptfmPitch from David
    infoStr=['PtfmPitch=' num2str(PtfmPitch_deg) 'deg'];
    
    % get data for plotting
    i_ptfmPitch=find(data_act.grid.PtfmPitch==PtfmPitch_deg);
    if ~isempty(i_ptfmPitch)
        % ptfmPitch_deg has been computed
        [x,y,z,isValid]=getGridData2(data_act,zName,xName,yName,'PtfmPitch',i_ptfmPitch);
    else
        % interpolation needed
        gridSmplPts_act={data_act.grid.BldPitchC data_act.grid.GenTq PtfmPitch_deg}; % {BldPitchC GenTq PtfmPitch}, see data_act.grid.dimNames
        [x,y]=ndgrid(data_act.grid.(xName),data_act.grid.(yName));
        z=squeeze(interpGridData(data_act,zName,gridSmplPts_act));
        isValid=squeeze(interpGridData(data_act,'isValid',gridSmplPts_act))==1;
    end

    % find maxima
    [zMax,xMax,yMax]=maxInterp2(z,x,y);
    [zMaxT,gridPtMax_arr]=maxInterpGridData(data_act,zName);

    % get ref
    GenTqRef_Nm=interp1(specs.schedule.v_mDs_arr,specs.schedule.GenTq_Nm_arr,data_act.windSpeed); 
    BldPitchCRef_deg=interp1(specs.schedule.v_mDs_arr,specs.schedule.BldPitch_deg_arr,data_act.windSpeed);

    % plot
    strTitle=[zName '@' infoStr '@wind=' num2str(data_act.windSpeed) 'm/s'];
    figure('NumberTitle','off','Name',strTitle)
    z(~isValid)=nan;
    x_iArr=1:size(x,1);
    y_iArr=1:size(y,2); % adapt if needed
    contourf(x(x_iArr,y_iArr),y(x_iArr,y_iArr),z(x_iArr,y_iArr),100,'LineStyle','none');
    hold on;
    plot(xMax,yMax,'rx','MarkerSize',10,'LineWidth',2);
%     plot(gridPtMax_arr(findselection(xName,data_act.grid.dimNames)),...
%          gridPtMax_arr(findselection(yName,data_act.grid.dimNames)),'or');
    xticks(x(:,1));
    yticks(y(1,:));
    xtickformat('%.2f')
    ytickformat('%.2f')
    xlabel(xName);
    ylabel(yName);
    title(zName);
    title(strTitle)
    hCB=colorbar;
    hCB.Label.String=zName;
    caxis([min(min(z(x_iArr,y_iArr))) max(max(z(x_iArr,y_iArr)))])
    grid on

    % mark and compare with reference value
    hline(gca,GenTqRef_Nm,'k-')
    vline(gca,BldPitchCRef_deg,'k-')
    GenPwrRef_W=interpGridData(data_act,'GenPwr',[BldPitchCRef_deg GenTqRef_Nm PtfmPitch_deg]);
    GenPwrImprvmnt_rel=(zMax-GenPwrRef_W)/GenPwrRef_W

end



%% plot different outputs over wind speed

close all
outName_arr={'GenPwr'                'TwrBsMyt'                           'GenSpeed'              'GenTq'                 'BldPitchC'          'PtfmPitch'                      'TSR'      'RotThrust'};
plotName_arr={'generator power (MW)' {'tower base' 'pitch moment  (MNm)'} 'generator speed (rpm)' 'generator torque (MNm)' 'blade pitch (deg)' {'platform' 'pitch angle (deg)'} 'TSR (-)'  'rotor thrust (MN)'};
plotFactor_arr=[1e-3            1e-3                                       1                       1e-6                    1                   1                                 1         1e-3];

for i_plot=[1 2  4 5 6 8] %1:numel(outName_arr)

    figure('Name',['Floating_' outName_arr{i_plot}],'NumberTitle','off');
    hold on
    plot([schedulePP.WindSpeed],...
         [schedulePP.(outName_arr{i_plot})]*plotFactor_arr(i_plot),...
         'LineStyle','-');

    ylabel(plotName_arr{i_plot})
    xlabel('wind speed (m/s)')
    grid on
    xlim(minmax([scheduleOpt.WindSpeed]))
    xticks(5:5:max(xlim))
    xtickangle(0);
    box on
    garyfyFigure
    
%     savePNG
%     saveLATEXPDF
end




%% compute AEP

% compute AEP
aepOpt=calcAEP([scheduleOpt.GenPwr],[scheduleOpt.WindSpeed]);
aepGT=calcAEP([scheduleGT.GenPwr],[scheduleGT.WindSpeed]);
aepPP=calcAEP([schedulePP.GenPwr],[schedulePP.WindSpeed]);

