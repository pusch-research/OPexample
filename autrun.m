OPexampleWorkDir=fileparts(mfilename('fullpath'));


% add OpenFAST-matlab toolbox
addpath(genpath(fullfile(OPexampleWorkDir, 'OpenFAST-matlab-toolbox')))
disp('> openfast matlab toolbox added..')

% add PMtools
run(fullfile(OPexampleWorkDir, 'PMtools','autrun.m'))

% add OPtools
run(fullfile(OPexampleWorkDir, 'OPtools','autrun.m'))

clear
