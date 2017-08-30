function dFluorPlotComboTraces( frameRate, actualF, baselineF, fRatio )
%Suitable to plot actual fluorescence, baseline fluorescence traces in
%single coordinate, and deltaF/F traces at the corresponding locations 
%on the other side.
%
% dFluorPlotComboTraces( frameRate, actualF, baselineF, fRatio )
%
%
%Input discription
%
% frameRate: Hz, frames per second
% actualF, baselineF, fRatio: arrays in which each row corresponds to one ROI



    [ numTraces, numFrames ] = size(actualF);
    t = 1/frameRate : 1/frameRate : numFrames/frameRate;

    figure('Color', 'w');

    for i = 1 : numTraces

        subplot(numTraces,2,i*2-1);
        plot1 = plot(t, baselineF(i,:), t, actualF(i,:), 'LineWidth', 2);
        set(plot1(1),'Color',[0.0784313753247261 0.168627455830574 0.549019634723663]);
        set(plot1(2),'Color',[0.847058832645416 0.160784319043159 0]);
        box off
        xlim([ t(1) t(end) ]);
        set(gca, 'YTick', [ ]);
        ylabel(['Cell ' num2str(i)], 'FontSize', 14)
        AddXunit(' s');

        subplot(numTraces,2,i*2);
        plot(t, fRatio(i,:)*100, 'LineWidth', 2, 'Color', [0 0.498039215803146 0]);
        box off
        xlim([ t(1) t(end) ]);
        ylabel('%');
        AddXunit(' s');

    end

    

end

