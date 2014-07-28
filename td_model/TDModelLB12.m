%% Simple Code for CSC/PR/MS TD Model with Actions.
% EAL 11/03/11
% Last updated: 04/25/12
% Does blocking and/or overshadowing. For L&B special issue.

%% Model Parameters
numms = 6;          %Number of microstimuli per stimulus
numstimuli = 1;     %Number of stimuli, including the reward/US
alpha = 0.05;       %Step-size
decay = 0.985;      %Memory Trace decay rate
gamma = 0.97;       %Discount factor
lambda = 0.95;      %Eligibility trace decay rate
sigma = 0.08;       %Microstimulus width
theta = 0.25;       %Action Threshold
upsilon = 0.9;      %Action Decay
presence = 0.2;     %Presence microstimulus level

%% Experiment Code:
% CS2 is used for Overshadowing and Blocking

%stimulustime(1) = 10;           %Conditioned Stimulus (CS) time
rewardt = 110;                   %Reward Time
stimulustime(1) = rewardt;       %Unconditioned Stimulus (US) time
%stimulustime(3) = 310;            %Conditioned Stimulus 2 (CS) time
%stimulusduration (1) = rewardt - stimulustime (1) -1;
stimulusduration (1) = 1;
%stimulusduration (3) = rewardt - stimulustime (3) -1;
numtrials = 200;                 %Last 3 trials are different probes    
numtimesteps = 300;
delay = true;                    %redundant
stimrep = 3;                     %Current code does not allow for mixed representations
%1 = CSC; 2 = Pres; 3 = MS
if stimrep == 1 %CSC             %Ensures that vectors are big enough to store data
    numms = numtimesteps-1;
end

%% Initialize Data Vectors:
x = zeros(1, (numms + 1)*numstimuli);       %Stimulus representation
w = zeros(1, (numms + 1)*numstimuli);       %Weight vector
e = zeros(1, (numms + 1)*numstimuli);       %Eligibility Traces
delta = zeros (numtrials, numtimesteps);    %TD Errors
value = zeros (numtrials, numtimesteps);    %Value functions
action = zeros (numtrials, numtimesteps);   %Response Levels
ms = zeros (numtimesteps, numms);           %Microstimulus levels
%maxaction = zeros (1, numtrials);           

%% Microstimulus Generation (NIPS2009 Version)
trace = 1;
for timestep = 1:numtimesteps
    for micro = 1:numms
        %ms (timestep, micro) = ((1/sqrt(2*pi))*exp(-((trace-(micro/numms))^2)/(2*(sigma^2))));
        ms (timestep, micro) = trace * ((1/sqrt(2*pi))*exp(-((trace-(micro/numms))^2)/(2*(sigma^2))));
    end
    trace = trace * decay;
end


%% Run Simuluation:
for trial = 1:numtrials
    oldvalue = 0;               % Reset on every trial
    %oldreward = 0;
    oldaction = 0;
    e = zeros(1, (numms + 1)*numstimuli);
    triallen = numtimesteps;
    
     %% Blocking Design (Cue Competition Set):
%     if trial < 201
%        stimulustime(3) = numtimesteps + 1; 
%     else
%        stimulustime(3) = 85;
%        stimulustime(1) = 85;
%        stimulusduration (1) = rewardt - stimulustime (1) -1;
%        stimulusduration (3) = rewardt - stimulustime (3) -1;
%     end
    
    rewardtime = rewardt;
    stimulustime(1) = rewardt;
    
    %% Intermingled Probes (Timing Set):
%     if mod (trial, 5) == 0
%         rewardtime = 500;
%         stimulustime(2) = rewardtime;
%         stimulusduration (1) = ((rewardt - stimulustime (1) -1) * 2);
%     else
%         stimulusduration (1) = rewardt - stimulustime (1) -1;
%     end
    
%     %% Probe Trials at End:
%     %Logic of this if statement is backward, but it works.
%     if trial == numtrials  % Last AB-
%          figure
%          rewardtime = numtimesteps + 1;
%          stimulustime(2) = numtimesteps + 1;  %No US
%          stimulustime(3) = 85;                %Bring Back CSB
%          alpha = 0;  
%     elseif trial == (numtrials - 1) %2nd to last: A-
%          figure
%          rewardtime = numtimesteps + 1;
%          stimulustime(2) = numtimesteps + 1;  %No US
%          stimulustime(1) = 210;                %Bring Back CSA
%          stimulustime(3) = numtimesteps + 1;  %No CSB
%          alpha = 0;  
%     elseif trial == (numtrials - 2) %3rd to last: B- 
%          figure
%          rewardtime = numtimesteps + 1;
%          stimulustime(2) = numtimesteps + 1;  %No US
%          stimulustime(1) = numtimesteps + 1;  %No CSA
%          alpha = 0;                           %No Learning  
%     end
    
    %% Shift probe (not used in current version of paper)
    %if trial == numtrials
    %    stimulustime (3) = stimulustime (3) - 20;
    %end
    %    stimulustime(2) = rewardtime;
    
    for timestep = 1:triallen
        if (timestep == rewardtime)
        %if (timestep > 59) && (timestep < 65) %for multi-step US (unused)
            reward = 1;
        else 
            reward = 0;
        end
        %% Create the Stimulus Representation
        %% CSC Code:
        if stimrep == 1
            x = zeros(1, (numms + 1)*numstimuli);
            for stimulus = 1:numstimuli
                if timestep >= stimulustime (stimulus) 
                    x ((stimulus-1)*numtimesteps + (timestep-stimulustime (stimulus))+1) = 1;                    
                end
            end
        elseif stimrep == 2
        %% Presence microstimulus:
            for stimulus = 1:numstimuli
                if and(timestep <= (stimulustime(stimulus) + stimulusduration(stimulus)), and((timestep >= stimulustime(stimulus)), delay))
                    x (numstimuli * numms + stimulus) = presence;
                else
                    x (numstimuli * numms + stimulus) = 0;
                end
            end
        elseif stimrep == 3
        %% MS Code:
            %Hybird Representation (unused)
%             for stimulus = 1:numstimuli
%                 if and(timestep <= (stimulustime(stimulus) + stimulusduration(stimulus)), and((timestep >= stimulustime(stimulus)), delay))
%                     x (numstimuli * numms + stimulus) = presence;
%                 else
%                     x (numstimuli * numms + stimulus) = 0;
%                 end
%             end
            for stimulus = 1:numstimuli
                if (timestep >= stimulustime(stimulus))
                    %noise = rand*0.02 -.01; (unused)
                    noise = 0;
                    fprintf('stimulus: %d, timestep: %d\n', stimulus, timestep);
                    fprintf('x update indices: %s\n', num2str((stimulus-1)*numms +1:(stimulus-1)*numms +numms));
                    fprintf('ms lookup row: %d\n\n', timestep-stimulustime(stimulus)+1);
                    
                    x ((stimulus-1)*numms +1:(stimulus-1)*numms +numms) = ms (timestep-stimulustime(stimulus)+1,:) + noise;     
                else
                    x ((stimulus-1)*numms +1:(stimulus-1)*numms +numms) = 0;
                    x (numstimuli * numms + stimulus) = 0; %not sure I get this... these are the +1 indices in the ms specification above... 1:18 for 6 ms per 3 stim. 19:21 updated here, zeroed out.
                end
            end
        end
        %% Value Calculation
           
        %fprintf('trial: %d, x: %s\n', trial, num2str(x));
        value (trial, timestep) = dot(x,w);
        
        %% Action Selection:
        action (trial, timestep) = upsilon * oldaction + max (0, oldvalue - theta);
        oldaction = action (trial, timestep);
        
        %% Learning Algorithm:
        delta (trial, timestep) = reward + (gamma * value (trial, timestep)) - oldvalue; %TD Learning
       
        w = w + (alpha * delta (trial, timestep) * e);
        e = x + (gamma * lambda * e);
        oldvalue = dot(x,w); % Or oldvalue = value. Difference in which weights are used.
        %oldreward = reward;
    end
    %% Real-time plotting:
%     maxaction(k) = max(max(action));
    plot(value (trial, :), 'r');
    hold on;
    plot(delta (trial, :), 'b');
    plot(action (trial, :), 'g');
    hold off;
    drawnow;
end

%% Extract Data for Cue Competition Set
probetrialab(stimrep, :) = action(numtrials, :);
probetriala(stimrep, :) = action(numtrials-1, :);
probetrialb(stimrep, :) = action(numtrials-2, :);

%% Extract Data for Timing Set
%probetrial(stimrep, :) = action(numtrials, :);
%probeaction = action(5:5:500,:);
%[maxaction(stimrep,:) peaktime(stimrep, :)] = max(probeaction');

%% Comparison Figure
 figure
% plot(action(numtrials-2,:), 'r');
 hold on;
% plot(action(numtrials-1,:), 'b');
 plot(action(numtrials,:), 'g');
 plot(value(numtrials,:), 'b');
 hold off;
