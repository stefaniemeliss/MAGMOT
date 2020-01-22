% % % % clear all but:
% % % clearvars -except scanner initialEyetrackerTest eyetracking practice debug dummymode dryrunscanner phantom secondKeyboard debug_extreme ...
% % %     preLearningRest postLearningRest postDoPostQuestionnaire version EXP_DIR ...
% % %     subject group orderNumber startBlock totalBlocks startTrial

% Clear the workspace
close all;
clearvars;
sca;

% Clear the workspace
close all;
clearvars;

% get parameter input (optional)llllllllllllt
prompt = {'scanner', 'initial eye tracker calibration', 'do eye tracking', 'practice', 'debug'};
defaults = {'1', '',  '1', '1', '0'};
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
dummymode = 0; % can be 0 or; for eye tracking

dryrunscanner = false; % participant is in the scanner using the button box, but there is no trigger as the scan is not running
phantom = false; % scanner is running, but without a participant in it
secondKeyboard = false; % a second keyBoard is attached that is not button box

debug_extreme = false;

% specify which bits of the protocol should be run 
% note: if practice, the scripts said these values to be false
% automatically
preLearningRest = false;
postLearningRest = true;
postDoPostQuestionnaire = true;

version = 'PsychToolBox_script';

try
    EXP_DIR = ['/Users/steffimeliss/Dropbox/Reading/PhD/Magictricks/fmri_study/' version '/'];
    cd(EXP_DIR)
catch
    EXP_DIR = ['/Users/stefaniemeliss/Dropbox/Reading/PhD/Magictricks/fmri_study/' version '/'];
    cd(EXP_DIR)
end



%----------------------------------------------------------------------
%                       Experimemtal setup
%----------------------------------------------------------------------

% add path to be able to use the functions saved in EXP dir
addpath(EXP_DIR);

% get parameter input (optional)
prompt = {'subject ID (two digits)', 'group (int, ext)', 'order number (1-25)', 'start block (set 2 or 3 if necessary)', 'number of blocks', 'start trial (1-12; adjust if necessary)'};
defaults = {'',  '', '', '1', '3', '1'};
answer = inputdlg(prompt, 'Experimental Setup Information', 1, defaults);
% now decode answer
[subject, group, orderNumber, startBlock, totalBlocks, startTrial] = deal(answer{:});
orderNumber = str2double(orderNumber);
startBlock = str2double(startBlock); % block to start with
totalBlocks = str2double(totalBlocks); % number of blocks in total
startTrial = str2double(startTrial); % which trial to start with - in case MATLAB breaks at a certain number
motivation = group;
fMRI = ['MAGMOT_' subject];
if strcmp(group, 'int')
    motivation = -1;
elseif strcmp(group, 'ext')
    motivation = 1;
else sca;
    disp('*** wrong input of group, start again ***');
    return
end

% get input for eye tracker
if eyetracking
    %     if ~dummymode
    if ~IsOctave
        commandwindow;
    else
        more off;
    end
    
    % STEP 1
    % Added a dialog box to set your own EDF file name before opening
    % experiment graphics. Make sure the entered EDF file name is 1 to 8
    % characters in length and only numbers or letters are allowed.
    if IsOctave
        edfFile = 'DEMO';
    else
        prompt = {'Enter tracker EDF file name (1 to 8 letters or numbers)'};
        dlg_title = 'Create EDF file';
        num_lines= 1;
        def     = {'MAGMOT'};
        answer  = inputdlg(prompt,dlg_title,num_lines,def);
        %edfFile= 'DEMO.EDF'
        edfFile = answer{1};
        fprintf('EDFFile: %s\n', edfFile );
    end
    %     end
end

% read in trial order file and select the one corresponding to orderNumer
allTrialOrders = csvread('TrialOrders_best_maxRep-4.csv', 1, 0); % this creates an array with each trialOrder as a column
currentTrialOrder = allTrialOrders(:,orderNumber); % selects the column of allTrialOrders corresponding to orderNumber
if practice
    orderNumber = 1;
    startBlock = 1; % block to start with
    totalBlocks = 1; % number of blocks in total
    startTrial = 1;
    
    preLearningRest = false;
    postLearningRest = false;
    postDoPostQuestionnaire = false;
    
    currentTrialOrder = [1 2];
    
end

% general inputs
betweenRatingFixation = 0.05;
timeoutCuriosity = 6-betweenRatingFixation;
timeoutAnswer = 6;
pre_stim_rest = 2;

if debug_extreme
    timeoutCuriosity = 3;
    timeoutAnswer = 3;
end

if debug
    timeoutCuriosity = 4;
    timeoutAnswer = 4;
end

%----------------------------------------------------------------------
%                       Screen setup
%----------------------------------------------------------------------

% Setup PTB with some default values
PsychDefaultSetup(2);

% Get the screen numbers
screens = Screen('Screens');

% Select the external screen if it is present, else revert to the native screen
screenNumber = max(screens);
% screenNumber = min(screens);

% Seed the random number generator.
rand('seed', sum(100 * clock));

% Define black, white and grey
white = WhiteIndex(screenNumber);
grey = white / 2;
black = BlackIndex(screenNumber);

% define scaling factors
if debug_extreme
    scalingFactor = 0.5;
elseif debug
    scalingFactor = 0.75;
    %     scalingFactor = 1;
else
    scalingFactor = 1;
end

% define screen size and font size
pixelsize = [0, 0, scalingFactor*1440, scalingFactor*900];

% Open the screen
Screen('Preference', 'SkipSyncTests', 1); %forgo syncTests
if debug || debug_extreme
    [window, windowRect] = PsychImaging('OpenWindow', screenNumber, black, pixelsize, 32, 2);
else
    [window, windowRect] = PsychImaging('OpenWindow', screenNumber, black, [], 32, 2);
end

% Set the blend funciton for the screen
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% Get the size of the on screen window in pixels
[screenXpixels, screenYpixels] = Screen('WindowSize', window);

% Get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(windowRect);

% Flip to clear
Screen('Flip', window);

% Query the frame duration
ifi = Screen('GetFlipInterval', window);
slack = Screen('GetFlipInterval', window)/2;

timingCorrection = slack;

% Set the text size and font
fontSizeBig = round(scalingFactor*50);
fontSizeSmall = round(scalingFactor*30);
Screen('TextSize', window, fontSizeBig);
Screen('TextFont', window, 'Courier');

% hide cursor
HideCursor();

% Query the maximum priority level
topPriorityLevel = MaxPriority(window);

% set up language settings
feature('DefaultCharacterSet','ISO-8859-1');
Screen('Preference','TextEncodingLocale', 'en_US.ISO8859-1');


%----------------------------------------------------------------------
%                       Timing Information
%----------------------------------------------------------------------

% Interstimulus interval time in seconds and frames
isiTimeSecs = 1;
isiTimeFrames = round(isiTimeSecs / ifi); %isiTimeFrames roughly equals 1 sec

% Numer of frames to wait before re-drawing
waitframes = 1;

%----------------------------------------------------------------------
%                       Keyboard information
%----------------------------------------------------------------------

% define trigger device and pptResponseDevice
% Determine id for the trigger response
if secondKeyboard
    deviceString = 'Magic Keyboard';
else
    deviceString='932'; %scanner device: both button box and trigger; old button box: deviceString='fORP Interface';
end
computerString='Apple Internal Keyboard / Trackpad'; %scanner device: both button box and trigger; old button box: deviceString='fORP Interface';
[id,name] = GetKeyboardIndices;% get a list of all devices connected

triggerDevice=0;
for i=1:length(name)%for each possible device
    if strcmp(name{i},computerString)%compare the name to the name you want
        experimenterKeyboard=id(i);%grab the correct id, and exit loop
        if length(name) == 1
            triggerDevice = experimenterKeyboard;
        end
    elseif scanner
        if strcmp(name{i},deviceString)%compare the name to the name you want
            triggerDevice=id(i);%grab the correct id, and exit loop
            break;
        end
        %     elseif dryrunscanner
        %         if strcmp(name{i},deviceString)%compare the name to the name you want
        %             triggerDevice=id(i);%grab the correct id, and exit loop
        %             break;
        %         end
    elseif ~scanner && ~secondKeyboard
        triggerDevice = experimenterKeyboard;
        break
    end
end


if triggerDevice==0%%error checking
    error('No device by that name was detected')
    sca;
end


% if taking responses from the button box then participantResponses should be scanner device (same as triggerDevice),
% otherwise can be -1 to take responses from keyboard devices
if phantom
    pptResponseDevice = -1;
else
    pptResponseDevice = triggerDevice;
end

% % use whatever attached keyboard in debug mode
% if debug || debug_extreme
%     pptResponseDevice = -1;
%     triggerDevice = -1;
%     experimenterKeyboard = -1;
% end

% define keyboard names
escape = KbName('ESCAPE'); % used during video presentation
quit = KbName('q'); % used during instructions
experimenterKey = KbName('s');

keysOfInterestButtonbox=zeros(1,256);
if dryrunscanner
    keysOfInterestButtonbox(KbName({'r', 'g', 'b', 'y'}))=1;
elseif practice
    keysOfInterestButtonbox(KbName({'h', 'j', 'k', 'l'}))=1;
elseif ~scanner || secondKeyboard
    keysOfInterestButtonbox(KbName({'h', 'j', 'k', 'l'}))=1; % I might have to change the response ratings accoringly!
else
    keysOfInterestButtonbox(KbName({'r', 'g', 'b', 'y'}))=1;
end

keysOfInterestTrigger=zeros(1,256);
keysOfInterestTrigger(KbName({'t'}))=1;

keysOfInterestEyeTracking=zeros(1,256);
keysOfInterestEyeTracking(KbName({'Return', 'ESCAPE', 'c', 'd', 'v'}))=1;

%----------------------------------------------------------------------
%                       Define instructions
%----------------------------------------------------------------------

greeting1 = 'Hello and welcome!\n\nThank you very much\nfor participating in the experiment.';
greeting2 = 'Do you feel comfortable?\n\nIf you would like to, you can have\na little wiggle to make\nyourself even more comfortable.\n\nFor the scanning, it is very important\nthat you do not move your head.\n\nSo please try to find a position\nthat is as convenient as possible.';
greeting_prac = 'Hello and welcome!\n\nThank you very much\nfor participating in the experiment.';

headscout = 'We are going to run\na short localizer sequence. \n\nPlease keep your head as still as possible\nand continue reading.';
reminder = 'To remind you:\n\nPlease keep your head as still as possible\n\nand do not cross your legs or your arms.';
wait = 'Please wait.\n\nSomeone will talk to you shortly.';
resting = 'With the next sequence,\nwe are measuring your brain activity at rest.\n\nThe scan will last for approximately 10 min.\n You will see a white screen.\n\nPlease keep your eyes open\nand simply look at the white screen.\nYou are allowed to blink as usual.\n\nPlease try to NOT think about anything at all.';
restingStart = 'The screen will turn white shortly.\n\nPlease keep your eyes open\nand try not to think about anything.';

question = 'Do you have any questions?\n\nPlease just ask.';
trigger = 'scanning is starting, waiting for trigger';
keyPress = '(press any key to continue)';
keyPressGreen = 'Plese press GREEN key (ring finger)\nto confirm that you read this statement.';

experimenterInput = 'EXPERIMENTER INPUT:\ncontinue or abort';

keyExit = '(press any key to EXIT)';

%-------
post1 = 'The task is done, GOOD JOB!\n\nThank you very much for completing it.';
quest1 = 'Thank you.\n\nWe are nearly done.';
quest2 = 'We will start the last scan now.\nThis is a structural image of your brain.\nThat means you do not have to do \nanything at all.\n\nIt will take approximately 6 minutes.';
quest3 = 'To prevent you from being too bored,\nwe have prepared a questionnaire.\n\nThis questionnaire is about\nyour opinion of the experiment.';
quest4 = 'Each question can be answered on a scale\nfrom 1 (definitely disagree) to 7 (definitely\nagree). Similar to the curiosity\nratings, you have to move the red number\nto the number reflecting your opinion.\n\nTo move it to the left,\nplease use your index finger (blue button).\n\nTo move it to the right,\nplease use your middle finger (yellow button).\n\nTo confirm your selection,\nplease use your pinkie (red button).';
t1Start = 'The structural scans\nand the questionnaire will start shortly';
endOfExperiment = 'The experiment is over now.\n\nMany thanks for your help!\n\nWe are waiting for the scan to finish\nand will move you out of the scanner shortly.';

%-------
answerRating = 'How many people (out of 100)\nare able to correctly figure out the solution?';
curiosityRating = 'Please rate how CURIOUS\nyou were while watching the trick.';
%-------
practice1 = 'We are  going to practice\nthe actual task.\n\n Do you feel ready?\n\nWe will show you the instructions.';
practice2 = 'In this experiment you will be presented\nwith a series of magic tricks.\nThe videos are without audio.\n\nYour task is to carefully watch the videos\nand try to figure out what has happened.';
practice3 = 'Before the start of each magic trick\nyou will see a fixation point.';
practice4 = 'Afterwards, you are asked to give\nan estimate of how many people (out of 100)\nare able to correctly\nfigure out the solution to the trick.\n\nPossible answers are the following:\n\n0 - 10 people\n11 - 20 people\n21 - 30 people\n31 or more people';
practice5 = 'In addition to that, we would like you\nto rate how curious you were while watching\nthe magic tricks on a scale \n\nfrom 1 (not curious at all) to 7 (very curious)';
practice6 = 'For each of the answers,\nyou have 6 seconds. ';
practice7 = 'To select the estimate you think\nis correct, you have to press\nthe corresponding button on the button box.\n\nYour INDEX finger is lies on the blue button\ncorresponding to the answer "0 to 10 people".\nYour MIDDLE finger is on the yellow button\ncorresponding to "11 to 20 people".\nYour RING finger lies on the green button\ncorresponding to "21 to 30 people"\nand your PINKIE lies on the red button\ncorresponding to "31 and more people".';
practice8 = 'For the curiosity rating, you have\nto move the red number \nto the number representing your curiosity.\n\nTo move it to the left,\nplease use your index finger (blue button).\n\nTo move it to the right,\nplease use your middle finger (yellow button).\n\nTo confirm your selection,\nplease use your pinkie (red button).';
practice9 = 'You will see both the answer and\nthe rating screen to show you\nhow it looks like.\n\nAfter you have indicated your answer and\nyour rating, the coloured ink will turn white\nand you simply wait for the task to continue.';
block_prac1 = 'You will see two magictricks\nduring the practice.\n\nIn the actual experiment,\nthere will be 36 magic tricks.';
block_prac2 = 'The practice starts NOW.\n\nYou are asked to estimate\nhow many people are able to correctly find\nthe solution to the magic trick.';
%-------
task0 = 'Great,\nthe eye tracker calibration was successful!';
task1 = 'Do you feel okay?\n\nNext, we are  going to do\nthe actual experiment.\n\n Do you feel ready?\n\nWe will show you the instructions again.\n\n It is going to be the same task\nwe practised earlier.';
fieldmap = 'While you are reading the instructions,\nwe are running the fieldmap.\nSo please keep your head as still as possible.\n\nThere will be a screen asking you\nto wait so that someone can talk to you.\nPlease do so once you get there.';
task2 = 'In this experiment you will be presented\nwith a series of magic tricks.\n\nYour task is to carefully watch the videos\nand try to figure out what has happened.';
task4 = 'Afterwards, you are asked to give\nan estimate of how many people (out of 100)\nare able to correctly\nfigure out the solution to the trick.\n\nPossible answers are the following:\n\n0 - 10 people\n11 - 20 people\n21 - 30 people\n31 or more people';
task5 = 'In addition to that, we would like you\nto rate how curious you were while watching\nthe magic tricks on a scale \n\nfrom 1 (not curious at all) to 7 (very curious)';
task6 = 'For each of the answers,\nyou have 6 seconds. ';
taskExt1 = 'We ask you to answer the question\n"how many people are able to find\nthe solution?" and you can get\n';
taskExt2 = 'To remind you:\n\nWe ask you to answer the question\n"how many people are able to find\nthe solution?" and you can get\n';
taskExt0 = ' an additional 50% bonus payment on top\nof your payment for both tasks (GBP 30.00)\nif you answer all questions correctly.\nThat meanseach correct answer\nis worth an additional GBP 0.80';
block1 = 'In total, you will see 36 magic tricks.\nThese will be presented in 3 blocks.\n\nThere will be two breaks in between\nso that you can rest and relax.\n\nPlease try not to move at all\nwhile you do the task.';
block2 = 'The experiment is ready to START.\n\nYou are asked to estimate\nhow many people are able to correctly find\nthe solution to the magic trick.';
taskStart = 'The fixation point will show up shortly.';
%-------
break1_1 = 'Thank you,\nthe first block of the task is finished.\n\nWELL DONE!';
break1_2 = 'Thank you,\nthe second block of the task is finished.\n\nWELL DONE!';
break2 = 'Take a break for as long as you need to.\n\nThe next part of the experiment\nwill start as soon as you are ready.\n\n The task is going to be\nthe same as in the previous block.';
if strcmp(group, 'ext') % extrinsic
    break3 = 'The experiment is ready to CONTINUE.\nYou are again asked to estimate\nhow many people are able to correctly\nfind the solution to the magic trick\nand you can get\n';
else
    break3 = 'The experiment is ready to CONTINUE.\n\nYou are again asked to estimate\nhow many people are able to correctly\nfind the solution to the magic trick.';
end
%-------
calibrationInstruction = 'We will run the calibration\nof the eye tracker now.\n\nPlease follow the dot on the screen\nwith your eye.';
calibrationScreen = 'Press RETURN (on either simulus or host computer)\nto toggle camera image\n\nPress ESC to output/record\n\nPress C to calibrate\n\nPress V to validate';
calibrationCheck = 'First, we need to do the calibration again.\n\nAfterwards, the task starts!';
driftCorrection = 'Before the start of the next block,\nwe need to check the eye tracker calibration.\nThere will be a fixation target in the\ncentre of the screen. Please look at it.';
%-------
q1 = 'It was fun to do the experiment.';
q2 = 'It was boring to do the experiment.';
q3 = 'It was enjoyable to do the experiment.';
q4 = 'I was totally absorbed in the experiment.';
q5 = 'I lost track of time.';
q6 = 'I concentrated on the experiment.';
q7 = 'The task was interesting.';
q8 = 'I liked the experiment.';
q9 = 'I found working on the task interesting.';
q10 = 'The experiment bored me.';
q11 = 'I found the experiment fairly dull.';
q12 = 'I got bored.';
q13 = 'I put a lot of effort into this.';
q14 = 'I did not try very hard\nto do well at this activity.';
q15 = 'I tried very hard on this activity.';
q16 = 'It was important to me to do well at this task.';
q17 = 'I did not put much energy into this.';
q18 = 'I did not feel nervous at all while doing this.';
q19 = 'I felt very tense while doing this activity.';
q20 = 'I was very relaxed in doing this experiment.';
q21 = 'I was anxious while working on this task.';
q22 = 'I felt pressured while doing this task.';
q23 = 'I tried to find out how many people\nwill be able to find the solution.';
q24 = 'I was able to see the magic tricks properly.';

% define instruction lists for loop %
% instrListStart = {greeting1, greeting2, headscout, reminder, wait, resting, question, reminder};
instrListStart = {greeting1, greeting2, headscout, resting, question, reminder};
instrListEyetracking = {calibrationInstruction, calibrationScreen};

% instructions depending on whether eye tracking can be done or not
if eyetracking
    instrListTask = {task0, task1, fieldmap, task2, practice4, practice5, practice6};
    instrListBlock1 = {wait, block1, block2, calibrationCheck, reminder};  % in between drift check eye tracker
    instrListBlock2 = {break1_1, wait, break2, question, break3, calibrationCheck, reminder};
    instrListBlock3 = {break1_2, wait, break2, question, break3, calibrationCheck, reminder};
else
    instrListTask = {task1, fieldmap, task2, practice4, practice5, practice6};
    instrListBlock1 = {wait, block1, block2, reminder};
    instrListBlock2 = {break1_1, wait, break2, question, break3,  reminder};
    instrListBlock3 = {break1_2, wait, break2, question, break3, reminder};
end

% different set of instructions for practice
if practice
    instrListTask = {greeting_prac, practice1, practice2, practice3, practice4, practice5, practice6, practice7, practice8, practice9};
    endOfExperiment = 'The practice is over now.\n\nDo you have any questions?';
    instrListBlock1 = {question, block_prac1, block_prac2}; %PRACTICE
end

instrListPost = {post1, wait, resting, question, reminder};
% instrListQuest = {quest1, quest2, quest3, quest4, question, reminder, t1Start};
instrListQuest = {quest1, quest2, reminder, quest3, quest4};
questionnaire = {q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, q14, ...
    q15, q16, q17, q18, q19, q20, q21, q22, q23, q24};
questionnaire = questionnaire(randperm(length(questionnaire)));


%----------------------------------------------------------------------
%                     Conditions and trials
%----------------------------------------------------------------------

% setting up a folder for each participant
SAVE_DIR = fullfile(pwd, 'behavioural_data', subject);
mkdir(SAVE_DIR);

% define stimulus folder depending on whether it is practice or not and clear potentially hidden files
if practice
    vidLoc= fullfile(EXP_DIR,'tricks', 'practice');
else
    vidLoc = fullfile(EXP_DIR,'tricks');
end
cd(vidLoc);
delete ._*

% create vid list
cd(EXP_DIR);
vidList = dir(fullfile(vidLoc, '*.mp4')); % find .mp4 files in the 'vid' directory
vidList = {vidList.name}; % we just need the filenames, this is our trial list
sort(vidList);% sort vidList ascending

% define number of trials in total and per block
numTrials = length(vidList);

if debug && ~practice
    numTrials = 6; %debug
elseif debug_extreme && ~practice
    numTrials = 3; %debug
end

trialList = numTrials/totalBlocks; % n in block
trialListBlock = trialList;

% jittering
if practice
    jitterVideo = [4 8];
    jitterRating = [6 4];
else
    filenameVideo = fullfile(EXP_DIR, 'videoJitter.tsv');
    fileID = fopen(filenameVideo);
    J = textscan(fileID, '%f');
    fclose(fileID);
    jitterVideo = J{1,1}';
    
    filenameRating = fullfile(EXP_DIR, 'ratingJitter.tsv');
    fileID = fopen(filenameRating);
    J = textscan(fileID, '%f');
    fclose(fileID);
    jitterRating = J{1,1}';
end

% post questionnaire
numQuestions = length(questionnaire);

% information to send marker: which events are important
eventListName = fullfile(EXP_DIR, 'whichTimings.tsv');
fileID = fopen(eventListName);
eventList = textscan(fileID, '%s %s %s %s %s %s %s %s %s %s %s');
fclose(fileID);
eventList = eventList{1,1}';

% information to send marker: what are the timestamps for each event?
timingsName = fullfile(EXP_DIR, 'allTimingsCell.tsv');
fileID = fopen(timingsName);
timingList = textscan(fileID, '%s %f %f %f %f %f %f %f %f %f %f %f');
fclose(fileID);
vidListTransposed = timingList{1,1};
timingListMat = cell2mat(timingList(2:size(timingList,2)));
timingListCell = num2cell(timingListMat);
timingListAll = [timingListCell, vidListTransposed];
timingListAllNamed = [timingListAll;eventList]; %timingListAllNamed{row, column}

% dynamically define number of tricks and events from the cell array created
numTricks = size(timingListAllNamed,1)-1;
nameTricksColumn = size(timingListAllNamed,2);
numEvents = size(timingListAllNamed,2)-1;
nameEventsRow = size(timingListAllNamed,1);
tricks = timingListAllNamed(1:numTricks,nameTricksColumn); % this refers to the column with the names of the magic tricks in it
events = timingListAllNamed(nameEventsRow,1:numEvents); % this refers to the row with the names of the events in it

%----------------------------------------------------------------------
%                     Make response matrixes
%----------------------------------------------------------------------

varsOfInterestTaskName = {'subject' 'fMRI' 'group' 'motivation' 'orderNumber' 'block' 'triggerTaskBlockRaw' 'startBlock' 'endBlock' 'stimID' ...
    'trial' 'whichVid' 'tTrialStart' 'tTrialEnd' 'durationTrial' 'fixationInitialDuration' 'pre_stim_rest'...
    'displayVidOnset' 'displayVidOffset' 'displayVidDuration'...
    'fixationPostVidOnset' 'fixationPostVidDuration' 'jitterVideo_trial' 'displayAnswerOnset' 'displayAnswerDuration' 'timeoutAnswer' 'responseAnswer' 'timestampAnswer' 'timestampAnswerWhite' 'rtAnswer' ...
    'fixationPostAnswerOnset' 'fixationPostAnswerDuration' 'betweenRatingFixation' 'displayCuriosityOnset' 'displayCuriosityDuration' 'timeoutCuriosity' 'responseCuriosity' 'timestampCuriosity' 'timestampCuriosityWhite' 'rtCuriosity'...
    'startValueCuriosity' 'clicksCuriosity' 'fixationPostCuriosityOnset' 'fixationPostCuriosityDuration' 'jitterRating_trial' };
varsOfInterestAllName = [varsOfInterestTaskName, events]; % combines elements of event list with varsOfInterestTask


varsOfInterestQuestName = {'subject' 'fMRI' 'group' 'motivation' 'question' 'text' 'startValueQuestion' 'clicksQuestion' 'answerQuestion' 'questionOnset' 'questionOffset' 'rtAnswer'};

% use this information to create
respMatTask = cell(numTrials, size(varsOfInterestAllName,2));
respMatQuest = cell(numQuestions, size(varsOfInterestQuestName,2));

%----------------------------------------------------------------------
%                       Experimental loop
%----------------------------------------------------------------------

if initialEyetrackerTest
    if eyetracking
        MAGMOT_EyeLinkInit
        fprintf('DEBUG: back to main script, after MAGMOT_EyeLinkInit\n');
    end
    
    Screen('TextSize', window, fontSizeBig);
    DrawFormattedText(window, 'Experimenter: press Q'    ,...
        'center', 'center', white);
    Screen('Flip', window);
    
    keyIsDown = 0;
    while ~keyIsDown %wating for any key press
        [keyIsDown,secs, keyCode] = KbCheck(experimenterKeyboard); %checking computer keyboard
        if keyCode(quit) == 1 % if q is pressed, experiment aborts
            Screen('Flip', window);
            Priority(0);
            KbStrokeWait;
            sca;
            initialEyetrackerTest = false;
            clearvars -except scanner initialEyetrackerTest eyetracking practice debug dummymode dryrunscanner phantom secondKeyboard debug_extreme ...
                preLearningRest postLearningRest postDoPostQuestionnaire version EXP_DIR ...
                subject group orderNumber startBlock totalBlocks startTrial
            
            fprintf('\n\n\n\n\n\initialEyetrackerTest is set to FALSE now\n\n\n\n\n\n');
            
        else
            keyIsDown = 0;
        end
    end
end


totalTrials = startTrial - 1;


%----------------------------------------------------------------------
%                       questionnaire and MPRAGE
%----------------------------------------------------------------------
if postDoPostQuestionnaire
    cd(EXP_DIR);
    block = 3;    
    
    % ---------------  instructions t1 and questionnaire  --------------- %
    for instr = 1:length(instrListQuest)
        Screen('TextSize', window, fontSizeBig);
        DrawFormattedText(window, instrListQuest{instr}    ,...
            'center', 'center', white);
        
        Screen('TextSize', window, fontSizeSmall);
        DrawFormattedText(window, keyPress,...
            'center', screenYpixels*0.9, white);
        Screen('Flip', window);
        
        KbQueueCreate(pptResponseDevice, keysOfInterestButtonbox);
        
        MAGMOT_waitForAnyKeyPress
        
    end
    
    
    % ---------------  questionnaire  --------------- %
    MAGMOT_postQuestionnaire
end
    
    %----------------------------------------------------------------------
    %                       save workspace and end experiment
    %----------------------------------------------------------------------
    
    datetimestamp = datestr(now, 'dd-mm-yyyy_HH:MM');
    
    filename = sprintf('%s_workspace_%s.mat',fMRI, datetimestamp);
    save(fullfile(SAVE_DIR,filename));
    
    
    % End of experiment screen. We clear the screen once they have made their
    % response
    Screen('TextSize', window, fontSizeBig);
    DrawFormattedText(window, endOfExperiment,...
        'center', 'center', white);
    Screen('TextSize', window, fontSizeSmall);
    DrawFormattedText(window, keyExit,...
        'center', screenYpixels*0.9, white);
    Screen('Flip', window);
    
    KbQueueCreate(pptResponseDevice, keysOfInterestButtonbox);
    KbQueueStart(pptResponseDevice);
    KbQueueWait(pptResponseDevice); % Wait until the 't' key signal is sent
    KbQueueRelease(pptResponseDevice);
    sca;
    
    if practice
        practice = false;
        fprintf('\n\n\n\n\n\npractice is set to FALSE now\n\n\n\n\n\n');
    end
    



