% this script concatenates fmri time series data AFTER pre-processing has been done
%

% necessary input: onset and duration of the event of intererest (the one
% that should stay in the newly created data set) in long format
clear all;
clc;


% define variables
TR = 2;
demean_data = {'no' 'yes'}; % this code includes a passage in which the average time series gets subtracted from time series
demean_data = {'no'}; % this code includes a passage in which the average time series gets subtracted from time series
ordered_data = {'no' 'yes'}; % should data be ordered?
ordered_data = {'yes'}; % should data be ordered?
iter = 0;
doSME = {'yes'};

% SME_outcome outcomes
SME_outcome = {'bothRemembered' 'bothForgotten'};

% define PATH: study drive study folder
PATH = '/storage/shared/research/cinn/2018/MAGMOT/';

% define BIDS dir which is used to define the subject cell
BIDS_DIR = [PATH 'MAGMOT_BIDS/'];

% define Nifti software path
software_dir = [PATH 'software/NIfTI_20140122/'];
addpath(software_dir);

% what is the root directory for the preprocessed data file
CONCAT_root = [PATH 'derivatives/magictrickwatching/concat/'];

% to define the subject list dynamically, we use the subject folder list in the concat directory
cd(BIDS_DIR);
delete .* % delete all invisible files
files = dir;
files = files(~ismember({files.name},{'.','..', 'scripts'})); % delete current and parent directory strings and the scripts dir
filenames = {files.name}; % from structure fils, only use the name
subdirs = filenames([files.isdir]); %only use those file names that correspond to folders
subjects = subdirs;
% subjects = {'sub-experimental016'}
% subjects_excl = { 'sub-experimental020' 'sub-experimental022' 'sub-experimental030'}; % not preprocessed yet
subjects_excl = { 'xxx' }; % not preprocessed yet

% for script creating purposes, we are only using the first 6 participants
% subjects = {'sub-control001' 'sub-control002' 'sub-control003' 'sub-experimental004' 'sub-experimental005' 'sub-experimental006'};
% subjects = {'sub-control001' 'sub-control002'};


% read in scan duration file (created with behavioural preprocessing r script)
scandurfile = fullfile(PATH,'derivatives', 'magictrickwatching', 'scripts', 'MAGMOT_informationAboutScanDuration.tsv');
fid = fopen(scandurfile);
scanDurations = textscan(fid, '%s %s %d %d %d %d', 'Delimiter','\t', 'HeaderLines', 1);
fclose(fid);

% loop for subjects
for s = 1:length(subjects)
    
    % there are four subjects with problematic data: 08, 09, 23, and 46
    % excluded for now
    if ~strcmp(subjects{s}, subjects_excl)
        SUBJ_ID = subjects{s};
        
        
        
        %%%% 1. step: combine all runs into one file %%%%
        
        % process input from scandurfile
        scanDurationsSubj = scanDurations{1,1}; % this cell has SUBJ_ID in it
        subj_match = strfind(scanDurations{1,1}, SUBJ_ID); % find rows of the cell matching the subject
        subj_index = find(not(cellfun('isempty', subj_match)));
        row_index = max(subj_index); % this is the row referring to last run of subject
        
        totalDurInTR =  scanDurations{1,6}(row_index); %this is the total duration of all three runs;
        fprintf('\n\n\n%s %s %d %s\n\n\n', SUBJ_ID, 'has total duration of', totalDurInTR, 'volumes');
        
        % define subject specific paths
        BIDS_DIR = [PATH 'MAGMOT_BIDS/' SUBJ_ID '/func/'];
        CONCAT_DIR = [CONCAT_root SUBJ_ID '/'];
        
        %%%% 1. step: concatenate the magictrick parts %%%%
        
        % open the nii file with the preprocessed data
        % note: after afni preprocessing, the final file already only
        % contains data of all three runs
        niifilename =[SUBJ_ID '_task-magictrickwatching_afniproc.nii.gz'];
        niifile = fullfile(CONCAT_DIR,niifilename);
        nii_subj = load_untouch_nii(niifile); % open nifti file
        img_subj = nii_subj.img; % use img data
        % feedback to console
        fprintf('%s %s %d %s\n', niifilename, 'has', size(img_subj,4), 'volumes');
        
        % get concat.tsv file containing information of magic tricks only
        cd(CONCAT_DIR);% go into concat folder
        concat_tsv = dir('*magictrickwatching*.tsv');
        eventfilename = {concat_tsv.name};
        eventfile = fullfile(CONCAT_DIR,eventfilename);
        
        % open events file
        fid = fopen(char(eventfile));
        eventsInSecs = textscan(fid, '%f %f %f %f %f %f %d %s %s %s %s %d %d', 'Delimiter','\t', 'HeaderLines', 1);
        % columns: vid (onset), mock (onset), duration_vid, duration_vid_withoutMock, duration_vid_withoutMock_postFixation, avgVidDur_MAGMOT,
        %            trial, stimFile, trialType_allConfRecognition, trialType_highConfRecognition, trialType_aboveAverageRecogntion, run, acq
        fclose(fid);
        
        % process input from eventfile
        stimIDs = eventsInSecs{1,8}; % get stimID
        [ordered_stimIDs, idx] = sort(stimIDs); % stimIDs start with two digits, so they can be sorted easily
        
        % create events array and convert seconds to volumes
        events(:,1) = round(eventsInSecs{1,2}/TR); % onsets: use second column as this reflects the end of the mock vid
        events(:,2) = ceil(round(eventsInSecs{1,6})/TR); % durations: use sixth column as this reflects average duration without mock vid
        
        for dd =  1:length(ordered_data) % create data in random and ordered way
            order = ordered_data{dd};
            
            % define onsets and duration
            onsets = events(:,1);
            durations = events(:,2);
            if strcmp(order,'yes')
                % create ordered arrays for onset and duration
                % idx rarranges the onset/duration array in that it reflects
                % when the magic trick first in the order of ordered_stimIDs began and how long it lasted
                onsets = onsets(idx);
                durations = durations(idx);
            end
            
            for ddd =  1:length(demean_data) % create data demeaned and not demeaned
                demean = demean_data{ddd};
                
                % cut and reorder the nifti file
                img_cut = [];
                count = 1; % to start with the first TR
                for ii = 1:size(onsets,1) % for all 36 magic tricks
                    ix = img_subj(:,:,:,onsets(ii):onsets(ii)+durations(ii)-1); % this is the time course for a defined magic trick
                    ixm = mean(ix,4); % this is the mean brain activation for a defined magic trick
                    ixdem = zeros(size(ix)); % 4D file (length of the magic trick) with ZEROS
                    if strcmp(demean,'yes')
                        for tt = 1:size(ix,4) % for each of the volumes/TRs
                            ixdem(:,:,:,tt) = double(ix(:,:,:,tt))-ixm; % substract the mean for the magic trick from each of the volumes for the magic trick
                        end
                    else
                        for tt = 1:size(ix,4) % for each of the volumes/TRs
                            ixdem(:,:,:,tt) = ix(:,:,:,tt);
                        end
                    end
                    img_cut(:,:,:,count:count+durations(ii)-1) = int16(ixdem); % add the time course for the defined magic trick to img_cut at the position 'count'
                    count = count+durations(ii); % update count so that the next snippet is added at the right position in the forth dimension
                end
                
                % save data
                nii_new=nii_subj;
                nii_new.hdr.dime.dim(5)=size(img_cut,4); % apply information from img to the nifti
                nii_new.img = img_cut;
                
                % define output file name
                if strcmp(order,'yes')
                    if strcmp(demean,'yes')
                        outputfilename = regexprep(niifilename,'.nii.gz','_demeaned_ordered'); % for this file, the magic tricks have not been reordered yet!
                    else
                        outputfilename = regexprep(niifilename,'.nii.gz','_ordered'); % for this file, the magic tricks have not been reordered yet!
                    end
                else
                    if strcmp(demean,'yes')
                        outputfilename = regexprep(niifilename,'.nii.gz','_demeaned_random'); % for this file, the magic tricks have not been reordered yet!
                    else
                        outputfilename = regexprep(niifilename,'.nii.gz','_random'); % for this file, the magic tricks have not been reordered yet!
                    end
                end
                outputfile = fullfile(CONCAT_DIR,outputfilename);
                
                % save new nii file
                save_untouch_nii(nii_new, outputfile);
                
                % gzip the new nii file and remove the unzipped file
                gzip([outputfile, '.nii']);
                delete([outputfile, '.nii']); 
                
                
                % feedback to console
                fprintf('%s %d\n','cutted from', size(img_subj,4))
                fprintf('%s %d\n','to', size(img_cut,4))
                fprintf('%s %s\n','saved as',outputfile)
                
                
            end
            
        end
        
        
        % clear some vars to not run into memory issues
        clear nii_new img_cut ix ixdem eventsInSecs events
        
    end
    
    %%%%%%%%%%%%%%%%%% TO CREATE THE FILES FOR THE SUBSEQUENT MEMORY EFFECT %%%%%%%%%%%%%%%%%%
    
    if strcmp(doSME,'yes')
        
        cd(CONCAT_root)
        
        allfiles_subj = dir(['**/*_' SUBJ_ID '*_SME_concat.tsv']);
        filenames_subj = {allfiles_subj.name}; % from structure fils, only use the name
        foldernames_subj = {allfiles_subj.folder}; % from structure fils, only use the name
        
        for f=1:length(allfiles_subj)
            
            file_subj = fullfile(foldernames_subj{f}, filenames_subj{f});
            
            % open events file
            fid = fopen(file_subj);
            eventsInSecs_SME = textscan(fid, '%f %f %f %f %s %s %s %s %f %s %s %s', 'Delimiter','\t', 'HeaderLines', 1);
            % columns: vid (onset), mock (onset), duration_vid, duration_vid_withoutMock,
            %         trialType_highAllRecognition, trialType_highConfRecognition, trialType_aboveAverageRecogntion,
            %         stimFile, avgVidDur_MAGMOT, behavParcel_allConf, behavParcel_highConf, behavParcel_aboveAvgConf
            fclose(fid);
            
            % process input from eventfile
            stimIDs = eventsInSecs_SME{1,7}; % get stimID
            stimIDs = eventsInSecs_SME{1,8}; % get stimID
            [ordered_stimIDs, idx] = sort(stimIDs); % stimIDs start with two digits, so they can be sorted easily
            
            % create events array and convert seconds to volumes
            events_SME(:,1) = round(eventsInSecs_SME{1,2}/TR); % onsets: use second column as this reflects the end of the mock vid
            events_SME(:,2) = ceil(round(eventsInSecs_SME{1,9})/TR); % durations: use ninth column as this reflects average duration without mock vid
            
            % define onsets and duration and order them using idx
            onsets_SME = events_SME(:,1);
            durations_SME = events_SME(:,2);
            onsets_SME = onsets_SME(idx);
            durations_SME = durations_SME(idx);
            
            % cut and reorder the nifti file
            img_cut = [];
            count = 1; % to start with the first TR
            for ii = 1:size(onsets_SME,1) % for all 36 magic tricks
                ix = img_subj(:,:,:,onsets_SME(ii):onsets_SME(ii)+durations_SME(ii)-1); % this is the time course for a defined magic trick
                ixm = mean(ix,4); % this is the mean brain activation for a defined magic trick
                ixdem = zeros(size(ix)); % 4D file (length of the magic trick) with ZEROS
                if strcmp(demean,'yes')
                    for tt = 1:size(ix,4) % for each of the volumes/TRs
                        ixdem(:,:,:,tt) = double(ix(:,:,:,tt))-ixm; % substract the mean for the magic trick from each of the volumes for the magic trick
                    end
                else
                    for tt = 1:size(ix,4) % for each of the volumes/TRs
                        ixdem(:,:,:,tt) = ix(:,:,:,tt);
                    end
                end
                img_cut(:,:,:,count:count+durations_SME(ii)-1) = int16(ixdem); % add the time course for the defined magic trick to img_cut at the position 'count'
                count = count+durations_SME(ii); % update count so that the next snippet is added at the right position in the forth dimension
            end
            
            % save data
            nii_new=nii_subj;
            nii_new.hdr.dime.dim(5)=size(img_cut,4); % apply information from img to the nifti
            nii_new.img = img_cut;
            
            
            % define output file name
            if strcmp(demean,'yes')
                outputfilename_SME = regexprep(filenames_subj{f},'_concat.tsv','_afniproc_demeaned'); % for this file, the magic tricks have not been reordered yet!
            else
                outputfilename_SME = regexprep(filenames_subj{f},'_concat.tsv','_afniproc');
            end
            outputfile_SME = fullfile(foldernames_subj{f},outputfilename_SME);
            
            % save new nii file
            save_untouch_nii(nii_new, outputfile_SME);
            fprintf('%s  %s %d %s\n', outputfilename_SME,'with', size(img_cut,4), 'volumes')
            
            % gzip the new nii file and remove the unzipped file
            gzip([outputfile_SME, '.nii']);
            delete([outputfile_SME, '.nii']);
            
            % clear some vars to not run into memory issues
            clear nii_new img_cut ix ixdem eventsInSecs_SME events_SME
            
            
        end
    end
end


% note: this code could be written far more efficiently IF rather than
% looping through the pairs, one would loop through the subject list and
% once a .nii is loaded, just do all computations with that file
