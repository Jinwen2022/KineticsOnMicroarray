function [ax1]= heatscatter_overlap_ver2(X, Y, proteinColor,numbins,lims, plot_colorbar,dilutionFactor,ax1,isSecondPlot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% heatscatter(X, Y, numbins,lims, ax,numToShow, markersize, marker, plot_colorbar, plot_lsf, xlab, ylab, title)
% mandatory:
%            X                  [x,1] array containing variable X
%            Y                  [y,1] array containing variable Y
% optional:
%            numbins            [double], default 50
%                                number if bins used for the
%                                heat3-calculation, thus the coloring
%            lims                [minimumX maximumX minimumY maximumY]
%                                when creatining the grid
%           isSecondPlot        [logical] decribing if this is the second
%                               plot for overlaping heatscatters 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% mandatory
    if ~exist('X','var') || isempty(X)
        error('Param X is mandatory! --> EXIT!');
    end
    if ~exist('Y','var') || isempty(Y)
        error('Param Y is mandatory! --> EXIT!');
    end
    
    %%%% optional
    if ~exist('numbins','var') || isempty(numbins)
        numbins = 30;
    else
        % force number, not char input
        numbins = str2double(numbins);
    end
    
    if ~exist('lims','var') || isempty(lims)
        [values, centers] = hist3([X(:) Y(:)], [numbins numbins]);
    else
        edges{1}=linspace(lims(1),lims(2),numbins);
        edges{2}=linspace(lims(3),lims(4),numbins);
       [values, centers] = hist3([X(:) Y(:)],'Edges',edges);
    end

    if ~exist('isSecondPlot','var') || isempty(isSecondPlot)
        isSecondPlot = false;
    end
    
    if ~exist('ax1','var') || isempty(ax1)
        figure; % Create a new figure
        ax1 = gca; % Get the current axes
    end
    numToShow = numel(X);
    markersize = 100*0.5587;
    marker = '.';
    centers_X = centers{1,1};
    centers_Y = centers{1,2};
    binsize_X = abs(centers_X(2) - centers_X(1)) / 2;
    binsize_Y = abs(centers_Y(2) - centers_Y(1)) / 2;
    bins_X = zeros(numbins, 2);
    bins_Y = zeros(numbins, 2);
    for i = 1:numbins
        bins_X(i, 1) = centers_X(i) - binsize_X;
        bins_X(i, 2) = centers_X(i) + binsize_X;
        bins_Y(i, 1) = centers_Y(i) - binsize_Y;
        bins_Y(i, 2) = centers_Y(i) + binsize_Y;
    end
    scatter_COL = zeros(length(X), 1);
    onepercent = round(length(X) / 100);
        
    for i = 1:length(X)
%         if (mod(i,onepercent) == 0)
%             fprintf('.');
%         end            
        last_lower_X = NaN;
        last_higher_X = NaN;
        id_X = NaN;
        c_X = X(i);
        last_lower_X = find(c_X >= bins_X(:,1),1,'last');
        if (isempty(last_lower_X))
            last_higher_X = find(c_X <= bins_X(:,2),1,'first');
        end
        if (~isnan(last_lower_X))
            id_X = last_lower_X;
        else
            if (~isnan(last_higher_X))
                id_X = last_higher_X;
            end
        end
        last_lower_Y = NaN;
        last_higher_Y = NaN;
        id_Y = NaN;
        c_Y = Y(i);
        last_lower_Y = find(c_Y >= bins_Y(:,1),1,'last');
        if (~isempty(last_lower_Y))
            last_higher_Y = find(c_Y <= bins_Y(:,2),1,'first');
        end
        if (~isnan(last_lower_Y))
            id_Y = last_lower_Y;
        else
            if (~isnan(last_higher_Y))
                id_Y = last_higher_Y;
            end
        end
        if ~isnan(id_X) && ~isnan(id_Y)
            scatter_COL(i) = values(id_X, id_Y);
        end
    
    end
    
    %Show density in the colorbar instead of number of operators
    scatter_COL=scatter_COL*(numToShow)/numel(X);
    
    standardAxesPosition = [0.1, 0.1, 0.7, 0.8]; % Adjust as needed

    % Determine if we are plotting the first or second set
    if ~isSecondPlot
        % This is the first plot
        set(ax1, 'Position', standardAxesPosition); % Set size to standard
    else
        % This is the second plot, overlay on the first
        ax2 = axes('Position', standardAxesPosition, 'Color', 'none'); % Ensure same position
        linkaxes([ax1, ax2],'xy'); % Link axes for zooming and panning synchronization
        ax1 = ax2; % Use the new axes for plotting
        hold on; % Retain existing plot
    end

    % Common scatter plot command for both
    scatter(ax1, X, Y, markersize, scatter_COL, 'Marker', marker, 'MarkerEdgeAlpha', 0.5, 'MarkerFaceAlpha', 0.5);

    % Set colormap and colorbar for the current axes
    numColors = 256;
    endColor = proteinColor;
%     dilutionFactor = 0.9;
    dilutedColor = endColor + dilutionFactor * (1 - endColor);
    startColor = dilutedColor;
    customColormap = [linspace(startColor(1), endColor(1), numColors)', ...
                      linspace(startColor(2), endColor(2), numColors)', ...
                      linspace(startColor(3), endColor(3), numColors)'];
    colormap(ax1, customColormap);

    if plot_colorbar
        cb = colorbar(ax1);
        if isSecondPlot
            cb.Position = [0.70, 0.5, 0.05, 0.3]; % Position of second colorbar
        else
            cb.Position = [0.60, 0.5, 0.05, 0.3]; % Position of first colorbar
        end
    end
   
    
    
end