
%% Combine Microstate and Behavioral Data Across Stability Analysis Runs
% This script loops over 5 clustering repetitions (each run of the stability analysis) and combines
% microstate statistics ("Microstats") with behavioral measures from the same participants.
%
% For each clustering run:
% 1. It merges microstate parameter data with HFS-related ratings (eg. intensity, unpleasantness).
% 2. It adds hypersensitivity area measures (length and width).
% 3. It links trial-level pinprick intensity and unpleasantness data, aligning with each participant's microstate results.
%
% Final outputs for each run are written to Excel files for later analysis.

%% Combine the Together data with Behavioural data
addpath(genpath('/Matlab Code'))
writepath = strcat('/StabilityAnalysis/');
run = ["Together1/","Together2/","Together3/","Together4/","Together5/"];

for thisrep = 1:5

    % For the clustering performed for the Experimental groups only
    locc = fullfile(writepath, run(thisrep));
    Microstats = readtable(strcat(locc,'Microstats.xlsx'));
    HFSint = readtable('HFSint.xls');
    HFSunp = readtable('HFSunp.xls');

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

    writetable(MicrostatsxHFS,strcat(locc,"MicrostatsxHFS.xlsx"))

    % For the clustering performed for the Experimental groups only
    clearvars -except writepath run locc
    Microstats = readtable(strcat(locc,'Microstats.xlsx'));
    Area = readtable('Area.xls');

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
    writetable(MicrostatsxArea,strcat(locc,"MicrostatsxArea.xlsx"))



    % For the clustering performed for the Experimental groups only
    clearvars -except writepath run locc
    Microstats = readtable(strcat(locc,'Microstats.xlsx'));
    Intensity = readtable('Intensity.xls');
    Unpleasantness = readtable('Unpleasantness.xls');
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

    writetable(MicrostatsxPPratings,strcat(locc,'MicrostatsxPPratings.xlsx'))

end



