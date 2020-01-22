% Clear the workspace
close all;
clearvars;
sca;

% get parameter input (optional)
prompt = {'scanner', 'initial eye tracker calibration', 'do eye tracking', 'practice', 'debug'};
defaults = {'1', '0',  '0', '', '0'};
title = 'Input: 1 (true) or 0 (false)';
dims = [1 50];
answer = inputdlg(prompt,title,dims,defaults);
% now decode answer
[scanner, initialEyetrackerTest, eyetracking, practice, debug] = deal(answer{:});
scanner = logical(str2num(scanner));
initialEyetrackerTest = logical(str2num(initialEyetrackerTest));
eyetracking = logical(str2num(eyetracking));
practice = logical(str2num(practice));
debug = logical(str2num(debug));

% ADDITIONAL parameters, most likely not to be changed
dummymode = 0; % can be 0 or; needed for eye tracking

dryrunscanner = false; % participant is in the scanner using the button box, but there is no trigger as the scan is not running
phantom = false; % scanner is running, but without a participant in it
secondKeyboard = false; % a second keyBoard is attached that is not button box

debug_extreme = false;

% specify which bits of the protocol should be run 
% note: if practice, the scripts said these values to be false
% automatically
preLearningRest = true;
doTask = true;
postLearningRest = true;
doPostQuestionnaire = true;

version = 'PsychToolBox_script';

try
    EXP_DIR = ['/Users/steffimeliss/Dropbox/Reading/PhD/Magictricks/fmri_study/' version '/'];
    cd(EXP_DIR)
catch
    EXP_DIR = ['/Users/stefaniemeliss/Dropbox/Reading/PhD/Magictricks/fmri_study/' version '/'];
    cd(EXP_DIR)
end

% then type: MAGMOT_PTB
MAGMOT_PTB

