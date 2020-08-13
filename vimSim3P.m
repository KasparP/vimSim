function [hF, axs]= vimSim3P(optsin)
%An abstract model of two-/three-photon imaging
%Simulates major physical constraints and allows optimization of imaging
%parameters

%written by Kaspar Podgorski, 2019-2020

opts.duration = 10;             tooltips.duration =     '(s) recording duration';
opts.reqFrameRate = 1e3;        tooltips.reqFrameRate = '(Hz) Minimum Frame rate';
opts.accessTime = 0;            tooltips.accessTime =   '(s) time cost per neuron to image multiple neurons. Assumes no access time cost for more pixels within a neuron.';
opts.dutyCycle = 1;             tooltips.dutyCycle =    'maximum fraction of time that the microscope records signals (generally 1 for Random Access imaging, <1 for line scanning)';
opts.power = 200;
opts.DoP = 1;                   tooltips.DoP =          'Degree of Parallelization; effective focus size, minimum 1';
opts.repRate = 1e8;             tooltips.repRate =      '(Hz) Number of laser pulses delivered to sample per second, before accounting for duty cycle and access time';
opts.nNeurons = 1;              tooltips.nNeurons =     'Number of neurons to image';
opts.indVar = fieldnames(opts); tooltips.indVar =   'independent variable to plot; leave empty for demo plots';
opts.doRandomAccess = false;    tooltips.doRandomAccess= 'only illuminate productive pixels';
opts.depth= 200;                tooltips.depth =        'imaging depth'; 
opts.figName = 'User-specified Parameters';
opts.do3P = false;              tooltips.do3P = 'Use 3P, rather than 2P, excitation';

if nargin %UPDATE WITH USER-SPECIFIED OPTIONS
    for field = fieldnames(optsin)'
         opts.(field{1}) = optsin.(field{1});
    end
    [hF, axs] = calculate(opts);
else
    hF = []; axs = [];
    optionsGUI(opts, tooltips);
end

end

function [hF1, axs] = calculate(opts)
opts = addHiddenParameters(opts);
[hF1, axs] = plotObj1(opts, opts.indVar);
set(hF1, 'name', opts.figName);
end

function opts = addHiddenParameters(opts)
%HIDDEN PARAMETERS
if opts.do3P
    opts.exP = 3;
    opts.exFac = 3e3;
    opts.bleachingP =4.5;            tooltips.bleachingP =   'power law for bleaching';
    opts.bleachingR =2e17;         tooltips.bleachingR =   'bleaching rate, arbitrary units';
    opts.attenuationLengthEx= 310;  tooltips.attenuationLengthEx = '(linear) attenuation length for excitation light'; %see e.g. Xu and Wise, 2013
    opts.damageP =4.5;                tooltips.damageP =      'power law for photodamage';
    opts.damageThresh= 1e-18;       tooltips.damageThresh = 'threshold for photodamage, arbitrary units';
else
    opts.exP = 2;
    opts.exFac = 1;
    opts.Sfac = 1e5;
    opts.bleachingP =3 ;            tooltips.bleachingP =   'power law for bleaching';
    opts.bleachingR =1e13;         tooltips.bleachingR =   'bleaching rate, arbitrary units';
    opts.damageP =3;                tooltips.damageP =      'power law for photodamage';
    opts.damageThresh= 1e-13;       tooltips.damageThresh = 'threshold for photodamage, arbitrary units';
    opts.attenuationLengthEx= 150;  tooltips.attenuationLengthEx = '(linear) attenuation length for excitation light'; %see e.g. Xu and Wise, 2013
end

opts.exSat = 1e-12;             tooltips.exSat =        'saturating excitation rate, arbitrary units';
opts.productiveDensity=0.02;    tooltips.productiveDensity = 'fraction of voxels that produce productive signal'; 
opts.nPix = 200;                tooltips.nPix =         'maximum number of pixels sampled per neuron';

opts.powLim = 150;              tooltips.powLim =       'heating threshold at point-sink limit, mW';
opts.powScale = 2.5e-4;         tooltips.powScale =     'heating threshold at large-area limit. mW/um^2 ';

%sanity Checks
opts.DoP = min(opts.DoP, opts.nPix);
opts.nNeurons = max(opts.nNeurons, 1);
opts.dutyCycle = max(0, min(1, opts.dutyCycle));
opts.accessTime = max(0, opts.accessTime);
end



function [hF, axs] = plotObj1(opts, indVar)

%find all fields with length>1
F = fieldnames(opts);
params = {};
np = 0;
for fix = 1:length(F)
    if length(opts.(F{fix}))>1 && isnumeric(opts.(F{fix}))
        np = np+1;
        pdata{np} = opts.(F{fix});
        params{np} = F{fix};
    end
end
if np==0
    error('no variable parameters were defined');
end
pgrid = cell(size(pdata));
[pgrid{:}] = ndgrid(pdata{:});
for pix = 1:length(pdata)
    opts.(params{pix}) = pgrid{pix}(:);
end

%CALCULATE OBJ FUNC 1
[S, effSat, damage, fracBleached, maxPower] = obj_signalRate(opts);
S = S./opts.nNeurons; %signal rate per neuron

%plots for each independent variable
if ~iscell(indVar)
    indVar = {indVar};
end
for indIx = 1:length(indVar)
    Sm=[]; Satm=[]; Damm=[]; Blm=[]; Pm=[];
    IV = indVar{indIx};
    if ~any(~strcmpi(params, IV))
        errordlg('Input a range for at least one of the dependent variables using matlab notation; examples: 1:100, logspace(0,3, 100)')
        error('Bad Input')
    end
    indVals = unique(opts.(IV));
    if length(indVals)<2
        errordlg('Input a range of values for the independent variable; examples: 1:100, logspace(0,3, 100)')
        error('Bad Input');
    end
    plotParams = params(~strcmpi(params, IV)); %we don't plot the independent variable
    hF = figure('units', 'normalized', 'pos', [0.2+0.1*indIx, 0.05, 0.3, 0.85]); axs = [];
    
    for ix = length(indVals):-1:1
        sel = opts.(IV)==indVals(ix);
        
        S_tmp = S(sel);
        [Sm(ix),I] = max(S_tmp);
        if ~isnan(Sm(ix))
            Sat_tmp = effSat(sel); Satm(ix) = Sat_tmp(I);
            Dam_tmp = damage(sel); Damm(ix) = Dam_tmp(I);
            Bl_tmp = fracBleached(sel); Blm(ix) = Bl_tmp(I);
            for pix = 1:length(plotParams)
                P_tmp = opts.(plotParams{pix})(sel);
                Pm(pix,ix) = P_tmp(I);
            end
        else
            Satm(ix) = nan;
            Damm(ix) = nan;
            Blm(ix) = nan;
            Pm(1:length(plotParams),ix) = nan;
        end
    end
    
    
    nplots = size(Pm,1)+6;
    colors = [0 0 0 ; (hsv(3)+1)./2 ; 0.5 0.5 0.5];
    xl = [min(indVals) max(indVals)];
    
    %plot independent variable
    axs(1) = subplot(nplots,1,1:3);
    plot(indVals, Sm, 'k', 'linewidth', 2);
    ylabel('Max Achievable Signal Rate Per Neuron');
    yl= round(max(Sm).*1.2, 2, 'significant');
    set(gca, 'ylim', [0 yl], 'xlim', xl);
    
    axs(2) = subplot(nplots,1,4); set(axs(2), 'UserData', 'sat'); %saturation
    plot(indVals, Satm, 'marker', '.', 'color', colors(2,:), 'linewidth',2);  ylabel('Saturation');
    set(gca, 'ylim', [0 1], 'xlim', xl);
    
    axs(3) = subplot(nplots,1,5); set(axs(3), 'UserData', 'bleach'); %bleaching
    plot(indVals, Blm, 'marker', '.', 'color', colors(3,:), 'linewidth',2); ylabel('Fraction Bleached');
    set(gca, 'ylim', [0 1], 'xlim', xl);
    
    axs(4) = subplot(nplots,1,6); set(axs(4), 'UserData', 'damage');%normalized damage
    plot(indVals, Damm, 'marker', '.', 'color', colors(4,:), 'linewidth',2); ylabel('Damage');
    set(gca, 'ylim', [0 1], 'xlim', xl);
    
    %plot auxiliary variables
    for pix = 1:length(plotParams)
        axs(pix+4) = subplot(nplots,1,pix+6); set(axs(pix+4), 'UserData', plotParams{pix});
        plot(indVals, Pm(pix,:), 'marker', '.', 'color', colors(5,:), 'linewidth',2); ylabel(plotParams{pix});
        switch plotParams{pix}
            case 'repRate'
                set(axs(end), 'yscale', 'log');
            case 'power'
                set(axs(end), 'ylim', [0 max(maxPower)*1.2]);
        end
    end
    xlabel(axs(end), IV, 'FontSize', 12);
    set(axs(end),'xlim', xl);
    
    switch IV
        case {'repRate', 'nNeurons'}
            set(axs, 'xscale', 'log')
    end
    linkaxes(axs, 'x')
end
end

function [S, exFrac, damageFrac, fracBleached, maxPower] = obj_signalRate(opts)
%OBJECTIVE FUNCTION
%computes the signal rate per neuron for given options

att  = exp(-opts.depth./opts.attenuationLengthEx); %attenuation factor
%compute the required area of the scan; this can be used to compute access time and heating constraints
if opts.doRandomAccess
    %adjust the effective rep rate according to duty cycle and access time to produce a long-run average:
    effRepRate = opts.repRate.*opts.dutyCycle.*max(0, (1-opts.accessTime.*opts.nNeurons.*opts.reqFrameRate));
    indRepRate = min(2e8, opts.repRate).*opts.dutyCycle.*max(0, (1-opts.accessTime.*opts.nNeurons.*opts.reqFrameRate)); %#independent measurements that can be made per second
    %adjust the number of pixels imaged per neuron, assuming even distribution of the available pulses
    opts.nPix = min(opts.nPix, opts.DoP .* indRepRate ./opts.reqFrameRate ./opts.nNeurons);
    bleachingPerVoxelPerPulse = opts.bleachingR.*(att.*opts.power./effRepRate./opts.DoP).^opts.bleachingP;
    bleachRate = (bleachingPerVoxelPerPulse./(1+bleachingPerVoxelPerPulse)) .* (opts.DoP ./ (opts.nNeurons .* opts.nPix));
else %raster scan
    effRepRate = opts.repRate.*opts.dutyCycle; %long-run average, not laser rate
    indRepRate = min(2e8, opts.repRate).*opts.dutyCycle;
    %adjust the number of pixels imaged per neuron, assuming even distribution of the available pulses
    opts.nPix = min(opts.nPix, opts.DoP .* indRepRate .*opts.productiveDensity./opts.reqFrameRate./opts.nNeurons);
    bleachingPerVoxelPerPulse = opts.bleachingR.*(att.*opts.power./effRepRate./opts.DoP).^opts.bleachingP;
    bleachRate = (bleachingPerVoxelPerPulse./(1+bleachingPerVoxelPerPulse)) .* opts.productiveDensity.* (opts.DoP ./ (opts.nNeurons .* opts.nPix));
end

recordingLengthInPulses = ceil(opts.duration .*effRepRate);
L = recordingLengthInPulses.*max(eps,bleachRate); %length of recording, in bleaching time constants
meanFracUnbleached = min(1, (1-exp(-L))./L); %mean brightness of the sample once bleaching is considered
fracBleached = 1-exp(-L); %

if opts.doRandomAccess
    ex = opts.exFac.*(att.*opts.power./effRepRate./opts.DoP).^opts.exP;  %excitation per focus per pulse, without saturation
    exSat = ex./(1 + ex./opts.exSat);        %incorporating saturation
    exFrac = exSat./opts.exSat;              %fraction of fluorophore excited with each laser pulse
    S = exSat .* opts.DoP .*effRepRate .*meanFracUnbleached;   %incorporating number of foci, number of pulses, and bleaching
    
    damagePerVoxelPerPulse =  (att.*opts.power./effRepRate./opts.DoP).^opts.damageP;   %mJ/pixel/pulse
    damage =  (opts.duration .* effRepRate .* opts.DoP ./ (opts.nNeurons .* opts.nPix) ) .* damagePerVoxelPerPulse; %average damage delivered to each pixel
    damageFrac = damage./opts.damageThresh;
else %raster
    ex = opts.exFac.*(att.*opts.power./effRepRate./opts.DoP).^opts.exP;  %excitation per focus per pulse, without saturation
    exSat = ex./(1 + ex./opts.exSat);        %incorporating saturation
    exFrac = exSat./opts.exSat;              %fraction of fluorophore excited with each laser pulse
    S = exSat .* opts.DoP .*effRepRate .*opts.productiveDensity.*meanFracUnbleached;   %incorporating number of foci, number of pulses, and bleaching
    
    damagePerVoxelPerPulse =  (att.*opts.power./effRepRate./opts.DoP).^opts.damageP;   %mJ/pixel/pulse
    damage =  (opts.duration .* effRepRate .* opts.DoP .*opts.productiveDensity ./ (opts.nNeurons .* opts.nPix) ) .* damagePerVoxelPerPulse; %average damage delivered to each pixel
    damageFrac = damage./opts.damageThresh;
end

%apply hard constraints
area = opts.nNeurons.*opts.nPix./opts.productiveDensity.*(0.4).^2; %in square microns
maxPower = sqrt((area.*opts.powScale).^2 + opts.powLim.^2); %simple model of maximum power vs area
S(indRepRate<opts.nNeurons.*opts.reqFrameRate) = nan;
S(opts.power>(maxPower+0.1)) = nan;
tooSlow = opts.nPix<1;
S(tooSlow) = nan;
S(damage>opts.damageThresh) = nan;

%Arbitrary Units
S = S .* 1e5;
end

function optsOut = optionsGUI(opts, tooltips)
%Podgorski Lab standard option GUI 2019
caller = dbstack;
if length(caller)>1
    caller = caller(2).name;
else
    caller = 'Unknown Function';
end
if nargin<2
    tooltips = [];
end

optsOut = opts;

[optNames, sortorder]= sort(fieldnames(opts)); %#ok<FLPST>
N = length(optNames);
titlesX = 5; titlesW = 80;
etX = titlesX+titlesW+20; etW = 200;
H = 15;
H0 =10;
HOK = 60;

handles.F = figure('Name', caller, 'pos', [600 600 titlesW+etW+40 N*(H+H0)+HOK+5], 'toolbar', 'none', 'menubar', 'none', 'resize', 'off', 'numbertitle', 'off');
handles.OK = uicontrol(...
    'Units','pixels',...
    'Parent',handles.F,...
    'Style','pushbutton',...
    'Position',[etX-50 5 etW 20],...
    'String','Demo Plots',...
    'Callback', @OKdemo);
handles.Demo = uicontrol(...
    'Units','pixels',...
    'Parent',handles.F,...
    'Style','pushbutton',...
    'Position',[etX 30 etW 20],...
    'String','Calculate',...
    'Callback', @OK);
handles.ET = [];

for n = 1:length(optNames)
    handles.titles(n) = uicontrol(...
        'Units','pixels',...
        'Parent',handles.F,...
        'Style','text',...
        'Position',[titlesX HOK+(n-1)*(H+H0)-4 titlesW H],...
        'String',optNames(n));
    
    handles.reset(n) = uicontrol(...
        'Units','pixels',...
        'Parent',handles.F,...
        'Style','pushbutton',...
        'Position',[etX+etW+5 HOK+(n-1)*(H+H0)-4 10 H],...
        'String','',...
        'Callback', @(varargin)(reset(n)));
    
    reset(n);
end
waitfor(handles.F);

    function parseET(src,n)
        type=class(opts.(optNames{n}));
        set(handles.titles(n), 'ForegroundColor', 'k')
        try
            switch type
                case 'logical'
                    optsOut.(optNames{n}) = eval(src.String{src.Value});
                case 'double'
                    optsOut.(optNames{n}) = eval(['[' src.String ']']);
                case 'cell'
                    optsOut.(optNames{n}) = eval(['''' src.String{src.Value} '''']);
                case 'char'
                    optsOut.(optNames{n}) = src.String;
            end
        catch ME
            %make the text error
            set(handles.titles(n), 'ForegroundColor', 'r')
        end
        
    end
    function OK(varargin)
        for k = 1:length(handles.ET)
            parseET(get(handles.ET(k)), k)
        end
        calculate(optsOut);
    end
    function OKdemo(varargin)
        for k = 1:length(handles.ET)
            parseET(get(handles.ET(k)), k)
        end
        demo(optsOut);
    end
    function reset(n)
        set(handles.titles(n), 'ForegroundColor', 'k')
        optsOut.(optNames{n}) = opts.(optNames{n});
        if length(handles.ET)>=n
            delete(handles.ET(n));
        end
        switch class(opts.(optNames{n}))
            case 'logical' %make a drop down menu
                handles.ET(n) = uicontrol(...
                    'Units','pixels',...
                    'Parent',handles.F,...
                    'Style','popupmenu',...
                    'String', {'false', 'true'},...
                    'Value', double(opts.(optNames{n}))+1,...
                    'Position',[etX HOK+(n-1)*(H+H0) etW H],...
                    'callback', @(src,evnt)(parseET(src, n)));
            case 'double'
                handles.ET(n) = uicontrol(...
                    'Units','pixels',...
                    'Parent',handles.F,...
                    'Style','edit',...
                    'String', num2str(opts.(optNames{n})),...
                    'Position',[etX HOK+(n-1)*(H+H0)-4 etW H],...
                    'callback', @(src,evnt)(parseET(src, n)));
            case 'char'
                handles.ET(n) = uicontrol(...
                    'Units','pixels',...
                    'Parent',handles.F,...
                    'Style','edit',...
                    'String', opts.(optNames{n}),...
                    'Position',[etX HOK+(n-1)*(H+H0)-4 etW H],...
                    'callback', @(src,evnt)(parseET(src, n)));
            case 'cell' %make a drop down
                handles.ET(n) = uicontrol(...
                    'Units','pixels',...
                    'Parent',handles.F,...
                    'Style','popupmenu',...
                    'String', opts.(optNames{n}),...
                    'Value', 1,...
                    'Position',[etX HOK+(n-1)*(H+H0) etW H],...
                    'callback', @(src,evnt)(parseET(src, n)));
            otherwise
                keyboard
        end
        if isfield(tooltips, optNames{n})
            set(handles.titles(n), 'TooltipString', tooltips.(optNames{n}));
            set(handles.ET(n), 'TooltipString', tooltips.(optNames{n}));
        end
    end
end