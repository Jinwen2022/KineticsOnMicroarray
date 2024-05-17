function heatscatter(X, Y, numbins,lims, ax,numToShow,markersize, marker, plot_colorbar, xlab, ylab, title)
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
%            numToShow          [double] Num of points to show in plot,
%                                don't plot the rest
%            ax                 handle to grafical axes for plotting
%            markersize         [double], default 10
%                                size of the marker used in the scatter-plot
%            marker             [char], default 'o'
%                                type of the marker used in the scatter-plot
%            plot_colorbar      [double], boolean 0/1, default 1
%                                set whether the colorbar should be plotted
%                                or not
%            xlab               [char], default ''
%                                lable for the x-axis
%            ylab               [char], default ''
%                                lable for the y-axis
%            title              [char], default ''
%                                title of the figure
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
    if ~exist('ax','var') || isempty(ax)
        figure
        ax = gca;
    end
    if ~exist('markersize','var') || isempty(markersize)
        markersize = 17;
    else
        % force number, not char input
        markersize = str2double(markersize);
    end
    if ~exist('marker','var') || isempty(marker)
        marker = '.';
    end
    if ~exist('plot_colorbar','var') || isempty(plot_colorbar)
        plot_colorbar = 1;
    end
    if ~exist('xlab','var') || isempty(xlab)
        xlab = '';
    end
    if ~exist('ylab','var') || isempty(ylab)
        ylab = '';
    end
    if ~exist('title','var') || isempty(title)
        title = '';
    end
    
    if ~exist('lims','var') || isempty(lims)
        [values, centers] = hist3([X(:) Y(:)], [numbins numbins]);
    else
        edges{1}=linspace(lims(1),lims(2),numbins);
        edges{2}=linspace(lims(3),lims(4),numbins);
       [values, centers] = hist3([X(:) Y(:)],'Edges',edges);
    end
    
    if ~exist('numToShow','var') || isempty(numToShow)
        numToShow = numel(X);
    end
    
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
    

%     figure();
%     ax=gcf;
    scatter(X(1:numToShow), Y(1:numToShow), markersize, scatter_COL(1:numToShow), marker);
    colormap(ax,cool)
    
    if (plot_colorbar)
        colorbar;
    end
    
    if (~isempty(xlab))
        xlabel(xlab);
    end
    if (~isempty(ylab))
        ylabel(ylab);
    end
    if (~isempty(title))
        title(title);
    end
    
end