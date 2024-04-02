function [specs,simpleModel]=specsUSFLOWTv3(FST_fileName,CpCqCt_fileName)

%% definitions

% file definitions
specs=struct();
specs.file.CpCqCt=CpCqCt_fileName;
specs.file.FST=FST_fileName;
specs.file.specs=mfilename("fullpath");

% meta data
specs.meta.string='USFLOWT-v3';
specs.meta.dateString=datestr(now,'YY-MM-DD_hh-mm');

% define operational parameters
specs.rGenPwr=10e3; %[kW] - units as in OpenFAST
specs.RotSpeed_min=2; % [rpm] - units as in OpenFAST
specs.WndSpeed_max=24; % cut off 
specs.WndSpeed_min=4; % cut in 


%% get all the data

% get data from FST file(s)
warning('off','FAST2Matlab:MoreHeadersThanCols');
warning('off','readFASTinFile:notSupported');
warning('off','readFASTinFile:noFile')
specs.FSTdata=readFASTinFiles(specs.file.FST,3);
warning('on','FAST2Matlab:MoreHeadersThanCols');
warning('on','readFASTinFile:notSupported');
warning('on','readFASTinFile:noFile')

% handle different versions OpenFAST input files
if isfield(specs.FSTdata,'ElastoFile')
    EDname='ElastoFile';
elseif isfield(specs.FSTdata,'EDFile')
    EDname='EDFile';
else
    error('not implemented.')
end


% build simpleModel
Jrot=156348016; % From ED.sum file
simpleModel=buildSimpleModel(...
    CpCqCt_fileName,...
    specs.rGenPwr,...
    specs.FSTdata.AeroFile.AirDens,...
    specs.FSTdata.(EDname).TipRad,...
    specs.FSTdata.(EDname).GBRatio,...
    Jrot,...
    {specs.WndSpeed_min specs.WndSpeed_max});

% copy easy access parameters
specs.R=specs.FSTdata.(EDname).TipRad;
specs.PreCone=specs.FSTdata.(EDname).PreCone_1_;
specs.ShftTilt=specs.FSTdata.(EDname).ShftTilt;
specs.GBRatio=specs.FSTdata.(EDname).GBRatio;
specs.Jrot=Jrot;
specs.Jgen=specs.FSTdata.(EDname).GenIner;
specs.Jtot=specs.Jrot+specs.Jgen*specs.GBRatio^2; % not fully correct
specs.NumBl=specs.FSTdata.(EDname).NumBl;
specs.rho=specs.FSTdata.AeroFile.AirDens;
specs.RtdRotSpeed=simpleModel.wr_rated_radDs*radDs2rpm;
specs.RtdGenSpeed=simpleModel.wg_rated_radDs*radDs2rpm;
specs.RtdGenTq=simpleModel.GenTq_rated_Nm; % [Nm] (input side)
specs.RtdWndSpeed=simpleModel.v_rated_mDs; 
specs.Cp_opt=simpleModel.Cp_opt;
specs.TSR_opt=simpleModel.TSR_opt;
specs.BldPitch_opt=simpleModel.BldPitch_opt_deg;
specs.schedule=simpleModel; % HACK: save whole simple model


%% old
% reference values from external source
% specs.ref.Cp_opt=0.4543;
% specs.ref.TSR_opt=9.2500;
% specs.ref.BldPitch_opt=1.3636;
% specs.ref.rGenTq=4.5403e+07; % [Nm]
% specs.ref.rWndSpeed=10.75; 
% specs.ref.TSR_opt=9.2500;
% specs.ref.BldPitch_opt=1.3636;
% specs.ref.



%% plot simple model

% surface
if nargout==0

    plotSimpleModel(simpleModel)

end

