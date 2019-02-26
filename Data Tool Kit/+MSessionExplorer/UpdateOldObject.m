function se = UpdateOldObject(se)
% 

if ischar(se) || iscellstr(se) || isstring(se)
    se = cellstr(se);
end

for i = 1 : numel(se)
    if iscellstr(se)
        % Load, update and save object
        s = load(se{i});
        structfun(@Update, s);
        save(se{i}, '-struct', 's');
    else
        % Update objet
        Update(se(i));
    end
end

end


function Update(se)

if ~isa(se, 'MSessionExplorer')
    return;
end

if isempty(se.epochInd) && se.numEpochs > 0
    disp('This object was generated before 2/26/2019');
    disp('Please refactor originalTrialInd to epochInd and numTrials to numEpochs in your code');
    se.epochInd = se.originalTrialInd;
else
    disp('This object is of the latest version');
end

end