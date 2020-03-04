%initialize matrices to hold all the x points, y points, and colors for plotting. The
%number of x and y points is equal to 3 times the the number of spikes for
%each cell.
all_xPoints = nan(sum(3*sum(cellfun(@length,all_near_laser_spikes_sorted_by_depth),2)),1);
all_yPoints = nan(sum(3*sum(cellfun(@length,all_near_laser_spikes_sorted_by_depth),2)),1);
all_colors = nan(sum(3*sum(cellfun(@length,all_near_laser_spikes_sorted_by_depth),2)),3);

figure()
index = 1;
for this_cell = 1:numCells-1
    % pull out the cell containing the spikes on each trial for a given cell 
    spikes = all_near_laser_spikes_sorted_by_depth(this_cell,:);
    % transpose the spikes (format required for plotSpikeRaster)
    spikes = cellfun(@transpose,spikes,'UniformOutput',false);
    numSpikes = sum(cellfun(@length,spikes));
    [xPoints, yPoints] = plotSpikeRaster(spikes,'PlotType','scatter');
    y_shift = size(all_near_laser_spikes_sorted_by_depth,2)*(numCells-this_cell);
    yPoints = yPoints + y_shift;
    this_cell_color = cmap(sorted_modulation_value(this_cell),:);
    
%     all_xPoints(index:index + 3*sum(cellfun(@length,spikes)) -1) = xPoints;
%     all_yPoints(index:index + 3*sum(cellfun(@length,spikes)) - 1) = yPoints;
%     all_colors(index:index + 3*sum(cellfun(@length,spikes)) - 1,:) =  this_cell_color;
    index = index + 3*sum(cellfun(@length,spikes));
    plot(xPoints,yPoints,'.k','MarkerSize',12,'color',this_cell_color)
    hold on
end



% plot(all_xPoints,all_yPoints,'LineStyle',':','LineWidth',0.35, 'color', all_colors)
% plot(all_xPoints,all_yPoints,'LineStyle','-','color','k','LineWidth',3)