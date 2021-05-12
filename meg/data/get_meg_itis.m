function get_meg_itis(tcfile)
load(tcfile);
ntrials=504;
iti_ideal=NaN(ntrials,1);
iti_onset=NaN(ntrials,1);
isi_ideal=NaN(ntrials,1);
feedback_onset=NaN(ntrials,1);
clock_onset=NaN(ntrials,1);
RTs=NaN(ntrials,1);
trial_offset=NaN(ntrials,1);
iti_offset=NaN(ntrials,1);
feedback_ideal=NaN(ntrials,1);
run=NaN(ntrials, 1);
trial=NaN(ntrials, 1);
rewFunc=cell(ntrials, 1);
emotion=cell(ntrials, 1);
score=NaN(ntrials, 1);
magnitude=NaN(ntrials, 1);
probability=NaN(ntrials, 1);
ev=NaN(ntrials,1);
isi_onset=NaN(ntrials, 1);
for i=1:ntrials
    if isempty(subject.timing(i).start)
        continue
    end
        
    if isfield(subject.timing(i), 'start')
        clock_onset(i) = subject.timing(i).start;
    end
    
    if isfield(subject.timing(i).face, 'onset')
        isi_onset(i) = subject.timing(i).face.onset;
    end
    
    %'iti' is actually the ISI between choice and feedback
    if isfield(subject.timing(i).ITI, 'expected')
        isi_ideal(i) = subject.timing(i).ITI.expected;
        feedback_onset(i) = subject.timing(i).ITI.onset;
    end
    
    if isfield(subject.timing(i).receipt, 'onset')
        iti_onset(i) = subject.timing(i).receipt.onset;      
    end
    
    if isfield(subject.timing(i).receipt, 'expected')
        feedback_ideal(i) = subject.timing(i).receipt.expected;
    end
    
    if isfield(subject.timing(i).ISI, 'expected')
        iti_ideal(i) = subject.timing(i).ISI.expected; %mislabeled
        iti_offset(i) = subject.timing(i).ISI.onset;
    end
       
    if isfield(subject.timing(i).end, 'onset')
        trial_offset(i) = subject.timing(i).end.onset;
    end
    
    rewFunc{i} = subject.order{i}{1};
    emotion{i} = subject.order{i}{12};
    run(i) = subject.order{i}{2};
    trial(i) = subject.order{i}{3};
    probability(i) = subject.order{i}{9};
    magnitude(i) = subject.order{i}{7};
    ev(i) = subject.order{i}{10};
    score(i) = subject.order{i}{8};
    RTs(i) = subject.order{i}{11}; %in ms
    
    %t_start is element 6 in subject.order, but encodes RT + 300ms in cumulative time (not relative to trial onset)
    %it is computed at the start of the scoreRxt function -- not esp. useful right now
    %isi_onset(i) = subject.order{i}{6};
    
    %the last trial in a block often has missingness in timings -- we can calculate it based on what's present
    if isnan(feedback_onset(i)), feedback_onset(i) = isi_onset(i) + .3; end %fixed duration
    if isnan(iti_offset(i)), iti_offset(i) = iti_onset(i) + subject.experiment{3}(i)/1000; end %timings hidden in subject.experiment
    if isnan(trial_offset(i)), trial_offset(i) = iti_onset(i) + subject.experiment{3}(i)/1000 + .02; end %usually about 20ms delay from iti offset to trial end
    if isnan(iti_ideal(i)), iti_ideal(i) = subject.experiment{3}(i)/1000; end %timings hidden in subject.experiment
    if isnan(isi_ideal(i)), isi_ideal(i) = .3; end %fixed duration
    if isnan(feedback_ideal(i)), feedback_ideal(i) = .85; end %fixed duration
end

subj_id = cell(ntrials, 1);
subj_id(:) = {subject.subj_id};

%isi_onset = isi_onset - clock_onset - RTs - .3; %hard code ISI ideal -- this is based on t_start. Not useful right now

xx=table(subj_id, run, trial, rewFunc, emotion, magnitude, probability, score,  RTs, clock_onset, isi_onset, feedback_onset, iti_onset, ...
    iti_offset, isi_ideal, feedback_ideal, iti_ideal, trial_offset, 'VariableNames', ...
    {'id', 'run', 'trial', 'rewFunc', 'emotion', 'magnitude', 'probability', 'score', 'rt', 'clock_onset', 'isi_onset', 'feedback_onset', ...
    'iti_onset', 'iti_offset', 'isi_ideal', 'feedback_ideal', 'iti_ideal', 'trial_offset'});

writetable(xx, strrep(tcfile, '.mat', '_alltimes.csv'));

end