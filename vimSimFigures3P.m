function vimSimFigures3P


%common parameters
opts.duration = 60;         tooltips.duration =     '(s) recording duration';
opts.reqFrameRate = 1e3;    tooltips.reqFrameRate = '(Hz) Minimum Frame rate';
opts.accessTime = 25e-6;    tooltips.accessTime =   '(s) time cost per neuron to image multiple neurons. Assumes no access time cost for more pixels within a neuron.';
opts.DoP = 1;               tooltips.DoP =          'Degree of Parallelization; effective focus size, minimum 1';
opts.nNeurons = 1;          tooltips.nNeurons =     'Number of neurons to image';
opts.damageP =3;            tooltips.damageP =      'power law for photodamage';
opts.depth= 500;

%% Panel 1: Performance and ideal rep rate vs imaging area

hF = [];
axs = {};

%raster, 2P
params = opts;
params.dutyCycle = 0.9;   
params.power = logspace(-1, log10(300), 300);
params.DoP = 1;    %degree of parallelization
params.repRate = logspace(5,10, 300);
params.nNeurons = 1:100;
params.indVar = 'nNeurons';
params.doRandomAccess = false;
params.do3P = false;
[hF(end+1), axs{end+1}] = vimSim3P(params); set(hF(end), 'name', 'raster 2p');
%raster, 3P
params.do3P = true;
[hF(end+1), axs{end+1}] = vimSim3P(params); set(hF(end), 'name', 'raster 3p');

%random access, 2P
params.DoP = 1;
params.repRate = logspace(5,10, 300);
params.doRandomAccess = true;
params.dutyCycle = 1;
params.do3P = false;
[hF(end+1), axs{end+1}] = vimSim3P(params); set(hF(end), 'name', 'random access 2P');
%random access, 3P
params.do3P = true;
[hF(end+1), axs{end+1}] = vimSim3P(params); set(hF(end), 'name', 'random access 3P');

%existing microscope
params.DoP = 3;
params.repRate= 8e7;
params.do3P = false;
params.nNeurons = 1:10;
[hF(end+1), axs{end+1}] = vimSim3P(params); set(hF(end), 'name', 'Existing microscope');

%merge relevant axes for figure
hF_merge = figure; %maximal signal rate
mergeAx1 = subplot(3,1,1:2); hold(mergeAx1, 'on'); ylabel('Signal Rate Per Neuron');
mergeAx2 = subplot(3,1,3); hold(mergeAx2, 'on'); ylabel('Ideal Rep Rate');
set(mergeAx1, 'xscale', 'log')
set(mergeAx2, 'xscale', 'log')

%set colors and styles
set(findobj(axs{1}, 'Type', 'line'), 'color', 'k');
set(findobj(axs{2}, 'Type', 'line'), 'color', 'b');

set(findobj(axs{3}, 'Type', 'line'), 'color', 'k', 'linestyle', ':', 'linewidth', 2);
set(findobj(axs{4}, 'Type', 'line'), 'color', 'b', 'linestyle', ':', 'linewidth', 2);

set(findobj(axs{5}, 'Type', 'line'), 'color', 'r');
perfThresh = get(findobj(axs{5}(1), 'Type', 'line'), 'ydata'); perfThresh = perfThresh(2);
plot(mergeAx1, [1 100], perfThresh.*ones(1,2), 'r', 'linestyle', '--', 'linewidth', 2);

for f_ix = [1 2 3 4]
    copyobj(get(axs{f_ix}(1), 'children'), mergeAx1); %signal
    copyobj(get(axs{f_ix}(6), 'children'), mergeAx2); %rep rate
end
copyobj(get(axs{5}(1), 'children'), mergeAx1); %extra

set(mergeAx2, 'yscale', 'log', 'xscale', 'log', 'ylim', [1e5 1e9])
set(mergeAx1, 'yscale', 'log');
xlabel(mergeAx2, '# of Neurons');
close(hF)


%% Panel 2 Performance and ideal rep rate vs depth
hF = [];
axs = {};

%raster
params = opts;
params.dutyCycle = 0.9;   
params.power = logspace(-1, log10(300), 300);
params.DoP = 1;    %degree of parallelization
params.repRate = logspace(4,10, 300);
params.nNeurons = 4;
params.depth = 0:20:1500;
params.indVar = 'depth';
params.doRandomAccess = false;
[hF(end+1), axs{end+1}] = vimSim3P(params); set(hF(end), 'name', 'raster 2p');
%raster, 3P
params.do3P = true;
[hF(end+1), axs{end+1}] = vimSim3P(params); set(hF(end), 'name', 'raster 3p');

%random access, 2P
params.DoP = 1;
params.repRate = logspace(4,10, 300);
params.doRandomAccess = true;
params.dutyCycle = 1;
params.do3P = false;
[hF(end+1), axs{end+1}] = vimSim3P(params); set(hF(end), 'name', 'random access 2P');
%random access, 3P
params.do3P = true;
[hF(end+1), axs{end+1}] = vimSim3P(params); set(hF(end), 'name', 'random access 3P');

%merge relevant axes for figure
hF_merge = figure; %maximal signal rate
mergeAx1 = subplot(3,1,1:2); hold(mergeAx1, 'on'); ylabel('Signal Rate Per Neuron');
mergeAx2 = subplot(3,1,3); hold(mergeAx2, 'on'); ylabel('Ideal Rep Rate');

%set colors and styles
set(findobj(axs{1}, 'Type', 'line'), 'color', 'k');
set(findobj(axs{2}, 'Type', 'line'), 'color', 'b');

set(findobj(axs{3}, 'Type', 'line'), 'color', 'k', 'linestyle', ':', 'linewidth', 2);
set(findobj(axs{4}, 'Type', 'line'), 'color', 'b', 'linestyle', ':', 'linewidth', 2);

for f_ix = [1 2 3 4]
    copyobj(get(axs{f_ix}(1), 'children'), mergeAx1); %signal
    if f_ix>2
        copyobj(get(axs{f_ix}(6), 'children'), mergeAx2); %rep rate
    end
end

plot(mergeAx1, [0 1500], perfThresh.*ones(1,2), 'r', 'linestyle', '--', 'linewidth', 2);

set(mergeAx2, 'yscale', 'log', 'ylim', [1e5 1e9])
set(mergeAx1, 'yscale', 'log', 'ylim', [1e-4 1]);
xlabel(mergeAx2, 'Depth');
close(hF)

end