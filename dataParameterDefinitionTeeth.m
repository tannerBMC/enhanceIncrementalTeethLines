%
% Definition of parameters and teeth data for extractCementumPatchesEnhanceIL program
%
% creates required output directories
%
% This program was used in 
%
% Tanner et al., Extended-field synchrotron microtomography for 
% non-destructive analysis of incremental lines in archeological human 
% teeth cementum, Proceedings of SPIE 11840 (2021) 1184019,
% DOI:10.1117/12.2595180
% 
% Christine Tanner, University of Basel, Biomaterials Science Center
%
% v1, 2024/04/15  cleaned up version

%% display parameters

showFigs = 1;      % 1: show main figure
showExtraFigs=0;   % 1: show additional figures

doSingle = 1;      % 1: do a single patch for a single slice
doScaleBars = 0;   % 1: do scale bar figures for manuscript
saveFigs = 1;      % 1: save main figures
saveExtraFigs = 0; % 1: save additional figures

saveData = 1;      % save results

% utility programs
path(path,['utils' filesep])

%% image data 

% tooth number T1-T4 in SPIE 2021 paper:
% (a) T1 & tooth4_hs04; hs07             2048 tif slices of size 11307x8507   idx 2
% (b) T2 & zahn_OKreC_probe1_hs02; hs04  2048 tif slices of size  8673x10773  idx 5
% (c) T3 & zahn_probe1_hs04              2048 tif slices of size 13990x11190  idx 4
% (d) T4 & zahn35_probe2_hs02            2048 tif slices of size 13066x9266   idx 3

sampleNamesV = {'doNotUse','tooth4','zahn35_probe2','zahn_probe1','zahn_OKreC_probe1'};
sampleIdx = 2;   % select sample index from 2, 3, 4, 5
sampleName = sampleNamesV{sampleIdx};
% select height step, 2, 4, 7 - not all available
if sampleIdx==2
    hs=7;   % 4 or 7 
elseif sampleIdx==3
    hs=2;   % 2
elseif sampleIdx==4
    hs=4;   % 4
elseif sampleIdx==5
    hs=2;   % 2 or 4
else
    disp([sampleName ' not available'])
end
    
hsStr = num2str(hs,'%02d');  % height step string
    
pixelDimensions=[0.65 0.65 0.65]; % image resolution [um]
maxSl=2048;                       % number of slices per height step

%% segmentation parameters

F=15;              % binning factor
distR = 30;        % ring size: 30 pixels from border in 15-binned image
lowW=40;           % lowW*F width in pixels for 15-binned image

%% patch and optimization parameters
     
% process projection in 3D with thickness thickness3D around given slice z
thickness3D = 35;         % >0: 3D thickness to consider [um]

% determine how many slices should be added on both sides
noSladd=round((thickness3D/pixelDimensions(3))/2);
noSl=2*noSladd+1;  

% patch parameters
pLen=50;   % patch length in pixels
years=70;  % defines size of cementum mask, assumes every year has size of 2um

% range of rotation to search, if not sufficient both increased by factor 1.5
maxRot=20.1;    % search in rotation range +-maxR [degree]
deltaRot=0.5;   % in steps of deltaRot [degree]

%% define teeth image data 

doErodeFill=0;  % default, 1: segmentation needs erosion/largestBody/fill

z=150;    % default start center slice for zV=z:2*noSl:maxSl-noSladd
if sampleIdx==2
    % 'tooth4'
    if hs==4
        z=55;   % start center slice 
        % select center slice and patch
        selZ=715;   % center slice for Fig. 6a and movie
        selP=17;    % patch for Fig. 6a and movie
    elseif hs==7
        z=330;   % start center slice 
        selZ=1650;  % for Fig. 6b and movie
        selP=27;    % 
    else
        disp([sampleStr ' ' hsStr ' not available'])
        return
    end
elseif sampleIdx==3
    % 'zahn35_probe2'
    if hs==2
        z=55;      % start center slice
        selZ=275;  % for Fig. 6f
        selP=14;   %
    else
        disp([sampleStr ' ' hsStr ' not available'])
        return
    end
elseif sampleIdx==4
    % 'zahn_probe1'
    if hs==4
        z=55;   % start slice
        selZ=275;   % for Fig. 6e
        selP=14;
    else
        disp([sampleStr ' ' hsStr ' not available'])
        return
    end
elseif sampleIdx==5
    % 'zahn_OKreC_probe1'
    doErodeFill=1;  % need erosion/largestBody/fill
    if hs==2
        z=880;   % start slice
        erodeR=2;
        selZ=1760;   % for Fig. 6c
        selP=35;
    elseif hs==4
        z=1320;   % start slice
        erodeR=10;
        selZ=1430;   % for Fig. 6d
        selP=8;
    else
        disp([sampleStr ' ' hsStr ' not available'])
        return
    end
else
    disp(['sample ' num2str(sampleIdx) ' not defined'])
    keyboard
end

subDirStr = [sampleName '_hs' hsStr];
if sampleIdx<4
    imageDir = ['dataDir1' filesep subDirStr filesep 'reco' filesep];
else
    imageDir = ['dataDir2' filesep subDirStr filesep 'reco' filesep];
end
outBaseDir = ['resultDir' filesep subDirStr filesep];

if ~exist(outBaseDir,'dir')
    mkdir(outBaseDir)
    disp(['created outBaseDir ' outBaseDir])
end

outDir = [outBaseDir filesep 'optRot_pLen' num2str(pLen) '_y' num2str(years) filesep];
if ~exist(outDir,'dir')
    mkdir(outDir)
    disp(['created outDir ' outDir])
end

% 3D processing
projStr = ['proj' num2str(thickness3D) '_'];

if doSingle
    % selected slice
    zV=selZ:selZ;
else
    % select which slices should be processed
    if (1==0)
        % process all slices
        zV=noSl:2*noSl:maxSl-noSladd;
    elseif (1==1)
        % process just single slice used for movie
        zV=selZ:selZ;
    else
        % process all slices starting from z
        zV=z:2*noSl:maxSl-noSladd;
    end
end
if length(zV)==1
    disp(['processing slab around slice ' num2str(zV(1))])
else
    disp(['processing slabs around ' num2str(length(zV)) ' slices in [' num2str(zV(1)) '-' num2str(zV(end)) ']'])
end
