
%% Combine microstate features with behavioral data for analysis
% This script merges microstate statistics (from EEG clustering) with behavioral outcomes
% from an experiment on HFS-induced hyperalgesia. The goal is to prepare a combined dataset
% for further statistical analysis, comparing EEG features with pain-related measures.

%% Combine the Together data with Behavioural data
% For the clustering performed for all the participants
writepath = 'path_to_all_microstate_results/';
Microstats = readtable(fullfile(writepath, 'Microstats.xlsx'));

datafolder = 'path_to_behavioral_microstate_data/';
HFSint = readtable(fullfile(datafolder, 'HFSint.xls'));
HFSunp = readtable(fullfile(datafolder, 'HFSunp.xls'));
Area = readtable(fullfile(datafolder, 'Area.xls'));
Intensity = readtable(fullfile(datafolder, 'Intensity.xls'));
Unpleasantness = readtable(fullfile(datafolder, 'Unpleasantness.xls'));

% HFSint and HFSunp
MicrostatsxHFS = Microstats(repelem(1:height(Microstats),5),:);
MicrostatsxHFS.HFSint = zeros(height(MicrostatsxHFS),1);
MicrostatsxHFS.HFSunp = zeros(height(MicrostatsxHFS),1);
MicrostatsxHFS.Train = zeros(height(MicrostatsxHFS),1);
MicrostatsxHFS.Fear = zeros(height(MicrostatsxHFS),1);
MicrostatsxHFS.Unpexp = zeros(height(MicrostatsxHFS),1);
MicrostatsxHFS.Intexp = zeros(height(MicrostatsxHFS),1);
participants = unique(Microstats.Subject);
for i = 1:length(participants)
    tempint = repmat(HFSint(HFSint.participant == participants(i),:),10,1);
    tempunp = repmat(HFSunp(HFSunp.participant == participants(i),:),10,1);
    participantidx = find(MicrostatsxHFS.Subject == participants(i));
    MicrostatsxHFS.HFSint(participantidx) = tempint.Ratings;
    MicrostatsxHFS.HFSunp(participantidx) = tempunp.Ratings;
    MicrostatsxHFS.Train(participantidx) = tempunp.Train;
    MicrostatsxHFS.Fear(participantidx) = tempunp.Fear;
    MicrostatsxHFS.Unpexp(participantidx) = tempunp.Unpexp;
    MicrostatsxHFS.Intexp(participantidx) = tempunp.Intexp;
end

writetable(MicrostatsxHFS,strcat(writepath,"MicrostatsxHFS.xlsx"))


% Area
MicrostatsxArea = Microstats;
MicrostatsxArea.Length = zeros(height(MicrostatsxArea),1);
MicrostatsxArea.Width = zeros(height(MicrostatsxArea),1);
participants = unique(Microstats.Subject);
for i = 1:length(participants)
    length = Area.area_sensitivity_length(Area.participant == participants(i));
    width = Area.area_sensitivity_width(Area.participant == participants(i));
    participantidx = find(MicrostatsxArea.Subject == participants(i));
    MicrostatsxArea.Length(participantidx) = length;
    MicrostatsxArea.Width(participantidx) = width;
end

writetable(MicrostatsxArea,strcat(writepath,"MicrostatsxArea.xlsx"))

% Pinprick Intensity and Unpleasantness
Intparticipant = zeros(height(Intensity),1);
Unpparticipant = zeros(height(Unpleasantness),1);
for i = 1:height(Intensity)
    subjectnoidx = isnumber(Intensity.ParticipantID{i});
    Intparticipant(i) = double(string(Intensity.ParticipantID{i}(subjectnoidx)));
    subjectnoidx = isnumber(Unpleasantness.ParticipantID{i});
    Unpparticipant(i) = double(string(Intensity.ParticipantID{i}(subjectnoidx)));
end
Intensity.ParticipantID = Intparticipant;
Unpleasantness.ParticipantID = Unpparticipant;

MicrostatsxPPratings = Microstats(repelem(1:height(Microstats),20),:);
MicrostatsxPPratings.Arm = strings(height(MicrostatsxPPratings),1);
MicrostatsxPPratings.IntensityRatings = zeros(height(MicrostatsxPPratings),1);
MicrostatsxPPratings.IntensityT0 = zeros(height(MicrostatsxPPratings),1);
MicrostatsxPPratings.UnpleasantnessRatings = zeros(height(MicrostatsxPPratings),1);
MicrostatsxPPratings.UnpleasantnessT0 = zeros(height(MicrostatsxPPratings),1);
MicrostatsxPPratings.Trial = zeros(height(MicrostatsxPPratings),1);
participants = unique(MicrostatsxPPratings.Subject);
for i = 1:length(participants)
    tempint = repmat(Intensity(Intensity.ParticipantID == participants(i),:),10,1);
    tempunp = repmat(Unpleasantness(Unpleasantness.ParticipantID == participants(i),:),10,1);
    participantidx = find(MicrostatsxPPratings.Subject == participants(i));
    MicrostatsxPPratings.Arm(participantidx) = string(tempint.Arm);
    MicrostatsxPPratings.IntensityRatings(participantidx) = tempint.IntensityRatings;
    MicrostatsxPPratings.IntensityT0(participantidx) = tempint.T0;
    MicrostatsxPPratings.UnpleasantnessRatings(participantidx) = tempunp.UnpleasantnessRatings;
    MicrostatsxPPratings.UnpleasantnessT0(participantidx) = tempunp.T0;
    MicrostatsxPPratings.Trial(participantidx) = tempint.Trial;
end
MicrostatsxPPratings.Arm(MicrostatsxPPratings.Arm == 'T1MFSint') = 'MFS';
MicrostatsxPPratings.Arm(MicrostatsxPPratings.Arm == 'T1CONTint') = 'Cont';

writetable(MicrostatsxPPratings,strcat(writepath,'MicrostatsxPPratings.xlsx'))


% For the clustering performed for the Control group only

writepath = 'path_to_microstate_results/';
Microstats = readtable(fullfile(writepath, 'Microstats.xlsx'));

datafolder = 'path_to_behavioral_data/';
HFSint = readtable(fullfile(datafolder, 'HFSint.xls'));
HFSunp = readtable(fullfile(datafolder, 'HFSunp.xls'));
Area = readtable(fullfile(datafolder, 'Area.xls'));
Intensity = readtable(fullfile(datafolder, 'Intensity.xls'));
Unpleasantness = readtable(fullfile(datafolder, 'Unpleasantness.xls'));

% HFSint and HFSunp
MicrostatsxHFS = Microstats(repelem(1:height(Microstats),5),:);
MicrostatsxHFS.HFSint = zeros(height(MicrostatsxHFS),1);
MicrostatsxHFS.HFSunp = zeros(height(MicrostatsxHFS),1);
MicrostatsxHFS.Train = zeros(height(MicrostatsxHFS),1);
MicrostatsxHFS.Fear = zeros(height(MicrostatsxHFS),1);
MicrostatsxHFS.Unpexp = zeros(height(MicrostatsxHFS),1);
MicrostatsxHFS.Intexp = zeros(height(MicrostatsxHFS),1);
participants = unique(Microstats.Subject);
for i = 1:length(participants)
    tempint = repmat(HFSint(HFSint.participant == participants(i),:),10,1);
    tempunp = repmat(HFSunp(HFSunp.participant == participants(i),:),10,1);
    participantidx = find(MicrostatsxHFS.Subject == participants(i));
    MicrostatsxHFS.HFSint(participantidx) = tempint.Ratings;
    MicrostatsxHFS.HFSunp(participantidx) = tempunp.Ratings;
    MicrostatsxHFS.Train(participantidx) = tempunp.Train;
    MicrostatsxHFS.Fear(participantidx) = tempunp.Fear;
    MicrostatsxHFS.Unpexp(participantidx) = tempunp.Unpexp;
    MicrostatsxHFS.Intexp(participantidx) = tempunp.Intexp;
end

writetable(MicrostatsxHFS,strcat(writepath,"MicrostatsxHFS.xlsx"))


% Area
MicrostatsxArea = Microstats;
MicrostatsxArea.Length = zeros(height(MicrostatsxArea),1);
MicrostatsxArea.Width = zeros(height(MicrostatsxArea),1);
participants = unique(Microstats.Subject);
for i = 1:length(participants)
    length = Area.area_sensitivity_length(Area.participant == participants(i));
    width = Area.area_sensitivity_width(Area.participant == participants(i));
    participantidx = find(MicrostatsxArea.Subject == participants(i));
    MicrostatsxArea.Length(participantidx) = length;
    MicrostatsxArea.Width(participantidx) = width;
end

writetable(MicrostatsxArea,strcat(writepath,"MicrostatsxArea.xlsx"))

% Pinprick Intensity and Unpleasantness
Intparticipant = zeros(height(Intensity),1);
Unpparticipant = zeros(height(Unpleasantness),1);
for i = 1:height(Intensity)
    subjectnoidx = isnumber(Intensity.ParticipantID{i});
    Intparticipant(i) = double(string(Intensity.ParticipantID{i}(subjectnoidx)));
    subjectnoidx = isnumber(Unpleasantness.ParticipantID{i});
    Unpparticipant(i) = double(string(Intensity.ParticipantID{i}(subjectnoidx)));
end
Intensity.ParticipantID = Intparticipant;
Unpleasantness.ParticipantID = Unpparticipant;

MicrostatsxPPratings = Microstats(repelem(1:height(Microstats),20),:);
MicrostatsxPPratings.Arm = strings(height(MicrostatsxPPratings),1);
MicrostatsxPPratings.IntensityRatings = zeros(height(MicrostatsxPPratings),1);
MicrostatsxPPratings.IntensityT0 = zeros(height(MicrostatsxPPratings),1);
MicrostatsxPPratings.UnpleasantnessRatings = zeros(height(MicrostatsxPPratings),1);
MicrostatsxPPratings.UnpleasantnessT0 = zeros(height(MicrostatsxPPratings),1);
MicrostatsxPPratings.Trial = zeros(height(MicrostatsxPPratings),1);
participants = unique(MicrostatsxPPratings.Subject);
for i = 1:length(participants)
    tempint = repmat(Intensity(Intensity.ParticipantID == participants(i),:),10,1);
    tempunp = repmat(Unpleasantness(Unpleasantness.ParticipantID == participants(i),:),10,1);
    participantidx = find(MicrostatsxPPratings.Subject == participants(i));
    MicrostatsxPPratings.Arm(participantidx) = string(tempint.Arm);
    MicrostatsxPPratings.IntensityRatings(participantidx) = tempint.IntensityRatings;
    MicrostatsxPPratings.IntensityT0(participantidx) = tempint.T0;
    MicrostatsxPPratings.UnpleasantnessRatings(participantidx) = tempunp.UnpleasantnessRatings;
    MicrostatsxPPratings.UnpleasantnessT0(participantidx) = tempunp.T0;
    MicrostatsxPPratings.Trial(participantidx) = tempint.Trial;
end
MicrostatsxPPratings.Arm(MicrostatsxPPratings.Arm == 'T1MFSint') = 'MFS';
MicrostatsxPPratings.Arm(MicrostatsxPPratings.Arm == 'T1CONTint') = 'Cont';

writetable(MicrostatsxPPratings,strcat(writepath,'MicrostatsxPPratings.xlsx'))












