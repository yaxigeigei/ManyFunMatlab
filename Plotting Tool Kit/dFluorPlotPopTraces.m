function [ shiftedFRatio ] = dFluorPlotPopTraces( frameRate, fRatio )
%Suitable to plot a population of deltaF/F traces in one panel with various colors
%
% [ shiftedFRatio ] = dFluorPlotPopTraces( frameRate, fRatio )
%
%
%Input discription
%
% frameRate: Hz, frames per second
% fRatio: an array in which each row corresponds to one ROI
%
%
%Output discription
%
% shiftedFRatio: the processed deltaF/F array used to plot the resulting figure


    [ numTraces, numFrames ] = size(fRatio);
    t = 1/frameRate : 1/frameRate : numFrames/frameRate;

    for i = 2 : numTraces
        shift = min(fRatio(i-1,:)) - 1.5 * (max(fRatio(i,:)) - min(fRatio(i,:)));
        fRatio(i,:) = fRatio(i,:) + shift;
    end
    
    figure('Color', 'w');
    plot(t, fRatio, 'LineWidth', 2);
    box off
    xlim([ t(1) t(end) ]);
    ylim([ min(min(fRatio)) max(max(fRatio)) ]);
    set(gca, 'YTick', [ ]);
    title('\Delta F / F', 'FontSize', 14);
%     AddXunit(' s');
    
    shiftedFRatio = fRatio;

    
    
end

