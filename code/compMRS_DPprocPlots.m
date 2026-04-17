clc
clear
close all

%% Load MAT file
load('CoMP-MRS-final.mat')   % loads variable: out

%% Settings
useRealPart = true;      % true = real(specs), false = abs(specs)
xRange = [0 4.5];
xTicks = 0:0.5:6.5;

% Requires boundedline on MATLAB path
% https://www.mathworks.com/matlabcentral/fileexchange/27485-boundedline-m

%% Check out size
[nRowsOut, nColsOut] = size(out);

if nRowsOut ~= 1 || nColsOut ~= 36
    warning('Expected out to be 1x36, but found %dx%d.', nRowsOut, nColsOut);
end

%% Find populated indices
validIdx = [];
for idx = 1:nColsOut
    if ~isempty(out{1,idx})
        validIdx(end+1) = idx; %#ok<AGROW>
    end
end

fprintf('Found %d populated out{1,index} entries.\n', numel(validIdx));
disp(validIdx)

%% Prepare tiled layout
nPlots = numel(validIdx);
nCols = ceil(sqrt(max(nPlots,1)));
nRows = ceil(nPlots / nCols);

figure('Color','w','Position',[100 100 1600 900]);
tiledlayout(nRows, nCols, 'TileSpacing','compact', 'Padding','compact');

plotSummary = struct();

%% Loop over each out{1,index}
for k = 1:nPlots
    idx = validIdx(k);
    dpLabel = sprintf('DP%02d', idx);
    thisBlock = out{1,idx};
    
    specMat = [];
    ppmRef = [];
    labels = strings(0,1);
    
    % Recursively collect spectra
    [specMat, ppmRef, labels] = collectSpectraFromBlock(thisBlock, useRealPart, specMat, ppmRef, labels);
    
    nexttile
    
    if isempty(specMat)
        title(sprintf('%s (no valid spectra)', dpLabel))
        axis off
        continue
    end
    
    meanSpec = mean(specMat, 1);
    stdSpec  = std(specMat, 0, 1);
    
    [hl, hp] = boundedline(ppmRef, meanSpec, stdSpec, 'b', ...
        'alpha', 'transparency', 0.3);
    
    xlim(xRange)
    xticks(xTicks)
    set(gca, 'XDir', 'reverse')
    xlabel('Frequency (ppm)')
    ylabel('Amplitude [a.u.]')
    title(sprintf('%s  (n = %d)', dpLabel, size(specMat,1)))
    
    ax = gca;
    ax.FontSize = 10;
    
    lgd = legend([hl hp], {'mean','\pm SD'}, 'Location','northwest');
    lgd.FontSize = 6;   % <<< smaller legend
    lgd.Box = 'off';
    
    box off
    
    plotSummary(k).index = idx;
    plotSummary(k).dpLabel = dpLabel;
    plotSummary(k).ppm = ppmRef;
    plotSummary(k).specMat = specMat;
    plotSummary(k).meanSpec = meanSpec;
    plotSummary(k).stdSpec = stdSpec;
    plotSummary(k).labels = labels;
end

saveas(gcf, 'boundedplots_each_out_index.png');
saveas(gcf, 'boundedplots_each_out_index.fig');

disp('Done: bounded plots created for each populated out{1,index}.');

%% ===== Local function =====
function [specMat, ppmRef, labels] = collectSpectraFromBlock(block, useRealPart, specMat, ppmRef, labels)

    if isempty(block)
        return
    end
    
    if iscell(block)
        for ii = 1:numel(block)
            [specMat, ppmRef, labels] = collectSpectraFromBlock(block{ii}, useRealPart, specMat, ppmRef, labels);
        end
        return
    end
    
    if isstruct(block)
        if isfield(block, 'ppm') && isfield(block, 'specs')
            try
                ppm = block.ppm(:).';
                
                if useRealPart
                    spec = real(block.specs(:)).';
                else
                    spec = abs(block.specs(:)).';
                end
                
                if isempty(ppmRef)
                    ppmRef = ppm;
                end
                
                if numel(ppm) == numel(ppmRef) && numel(spec) == numel(ppmRef)
                    specMat = [specMat; spec]; %#ok<AGROW>
                    
                    if isfield(block, 'filepath') && ~isempty(block.filepath)
                        try
                            labels(end+1,1) = string(block.filepath); %#ok<AGROW>
                        catch
                            labels(end+1,1) = "spectrum_" + string(size(specMat,1)); %#ok<AGROW>
                        end
                    else
                        labels(end+1,1) = "spectrum_" + string(size(specMat,1)); %#ok<AGROW>
                    end
                end
            catch
                % skip malformed entries
            end
        end
        return
    end
end