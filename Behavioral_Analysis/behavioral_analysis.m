clear all
close all
clc

%% Behavioral Data analysis

% 02/15/2023
% Daniel Sazhin

% This script takes a participant's responses and collapses them into a
% dollar amount for analysis.

% Update the subjects list with up-to-date subjects

maindir = pwd; % set on computer doing the analysis


mydir = pwd;
idcs = strfind(mydir,'\');
maindir = mydir(1:idcs(end)-1);
outputdir = '\Behavioral_Analysis\Behavioral_Results\';
output_folder = [maindir outputdir];

subjects = [2002];
%subjects_debug = [1240,1248]; % These have different column numbers, probably because RAs used psychopy2.

%% Extract data

% I need to parse the data into the separate games.

for ii= 1:length(subjects)
    subj = subjects(ii);
    
    UG_R_accept_2 = [];
    UG_R_reject_2 = [];
    DG_P_2 = [];
    UG_P_2 = [];
    UG_R_accept = [];
    UG_R_reject = [];
    DG_P = [];
    UG_P = [];

    try

r1 = 1; % Run 1
r2 = 1; % Run 2

% Select 1 if you want to use the run. Zero if not.

%% Run 1

if r1 > 0
    
    for r = 0
        
        % sub-101_task-ultimatum_run-0_raw.csv sub-102_task-ultimatum_run-1_raw.csv
        fname = fullfile(maindir,'logs',num2str(subj),sprintf('sub-%04d_task-ultimatum_run-%d_practice_raw.csv',subj,r)); % Psychopy taken out from Logs to make work for now.
        if exist(fname,'file')
            fid = fopen(fname,'r');
        else
            fprintf('sub-%d -- Let''s Make a Deal Game, Run %d: No data found.\n', subj, r+1)
            continue;
        end
        C = textscan(fid,repmat('%f',1,23),'Delimiter',',','HeaderLines',1,'EmptyValue', NaN);
        fclose(fid);
        
        
        % "Feedback" is the offer value (out of $20)
        
        Trial = C{1};
        Block = C{3};
        Endowment = C{4};
        response = C{17};
        response = round(response);
        L_Option = C{7};
        R_Option = C{8};
        
    end
    
    
    % Distribute games per responses
    
    %  # 4 = faculty camp 3 = students camp 2 = Penn faculty 1 = Penn student
        
    
     for t = 1:length(Endowment)
     
        % Block 2 = DG prop
        
        if Block(t) == 4 % If DG Proposer
            if response(t) == 2 % If Left button selected
                if round(L_Option(t) > R_Option(t)) % Is Left button more?
                    update = [Trial(t), Endowment(t), L_Option(t), 1, L_Option(t)/Endowment(t), R_Option(t)/Endowment(t)];  % Trial number, Endowment, Choice, More, M_Option_Proportion, L_Option_Proportion
                    DG_FC = [DG_FC; update];
                else % Therefore it is less.
                    update = [Trial(t), Endowment(t), L_Option(t), 0, R_Option(t)/Endowment(t), L_Option(t)/Endowment(t)];
                    DG_FC = [DG_FC; update];
                end
            end
        end
        
        if Block(t) == 4
            if response(t) == 3 % Right option
                if round(R_Option(t) > L_Option(t)) % Is right option more?
                    update = [Trial(t), Endowment(t), R_Option(t), 1, R_Option(t)/Endowment(t), L_Option(t)/Endowment(t)];  % Trial number, Endowment, Choice, More, M_Option_Proportion, L_Option_Proportion
                    DG_FC = [DG_FC; update];
                else
                    update = [Trial(t), Endowment(t), R_Option(t), 0, L_Option(t)/Endowment(t), R_Option(t)/Endowment(t)];
                    DG_FC = [DG_FC; update];
                end
            end
        end
        
        
         if Block(t) == 3 % If DG Proposer
            if response(t) == 2 % If Left button selected
                if round(L_Option(t) > R_Option(t)) % Is Left button more?
                    update = [Trial(t), Endowment(t), L_Option(t), 1, L_Option(t)/Endowment(t), R_Option(t)/Endowment(t)];  % Trial number, Endowment, Choice, More, M_Option_Proportion, L_Option_Proportion
                    DG_SC = [DG_SC; update];
                else % Therefore it is less.
                    update = [Trial(t), Endowment(t), L_Option(t), 0, R_Option(t)/Endowment(t), L_Option(t)/Endowment(t)];
                    DG_SC = [DG_SC; update];
                end
            end
        end
        
        if Block(t) == 3
            if response(t) == 3 % Right option
                if round(R_Option(t) > L_Option(t)) % Is right option more?
                    update = [Trial(t), Endowment(t), R_Option(t), 1, R_Option(t)/Endowment(t), L_Option(t)/Endowment(t)];  % Trial number, Endowment, Choice, More, M_Option_Proportion, L_Option_Proportion
                    DG_SC = [DG_SC; update];
                else
                    update = [Trial(t), Endowment(t), R_Option(t), 0, L_Option(t)/Endowment(t), R_Option(t)/Endowment(t)];
                    DG_SC = [DG_SC; update];
                end
            end
        end
        
         if Block(t) == 2 % If DG Proposer
            if response(t) == 2 % If Left button selected
                if round(L_Option(t) > R_Option(t)) % Is Left button more?
                    update = [Trial(t), Endowment(t), L_Option(t), 1, L_Option(t)/Endowment(t), R_Option(t)/Endowment(t)];  % Trial number, Endowment, Choice, More, M_Option_Proportion, L_Option_Proportion
                    DG_FP = [DG_FP; update];
                else % Therefore it is less.
                    update = [Trial(t), Endowment(t), L_Option(t), 0, R_Option(t)/Endowment(t), L_Option(t)/Endowment(t)];
                    DG_FP = [DG_FP; update];
                end
            end
        end
        
        if Block(t) == 2
            if response(t) == 3 % Right option
                if round(R_Option(t) > L_Option(t)) % Is right option more?
                    update = [Trial(t), Endowment(t), R_Option(t), 1, R_Option(t)/Endowment(t), L_Option(t)/Endowment(t)];  % Trial number, Endowment, Choice, More, M_Option_Proportion, L_Option_Proportion
                    DG_FP = [DG_FP; update];
                else
                    update = [Trial(t), Endowment(t), R_Option(t), 0, L_Option(t)/Endowment(t), R_Option(t)/Endowment(t)];
                    DG_FP = [DG_FP; update];
                end
            end
        end
        
         if Block(t) == 1 % If DG Proposer
            if response(t) == 2 % If Left button selected
                if round(L_Option(t) > R_Option(t)) % Is Left button more?
                    update = [Trial(t), Endowment(t), L_Option(t), 1, L_Option(t)/Endowment(t), R_Option(t)/Endowment(t)];  % Trial number, Endowment, Choice, More, M_Option_Proportion, L_Option_Proportion
                    DG_SP = [DG_SP; update];
                else % Therefore it is less.
                    update = [Trial(t), Endowment(t), L_Option(t), 0, R_Option(t)/Endowment(t), L_Option(t)/Endowment(t)];
                    DG_SP = [DG_SP; update];
                end
            end
        end
        
        if Block(t) == 1
            if response(t) == 3 % Right option
                if round(R_Option(t) > L_Option(t)) % Is right option more?
                    update = [Trial(t), Endowment(t), R_Option(t), 1, R_Option(t)/Endowment(t), L_Option(t)/Endowment(t)];  % Trial number, Endowment, Choice, More, M_Option_Proportion, L_Option_Proportion
                    DG_SP = [DG_SP; update];
                else
                    update = [Trial(t), Endowment(t), R_Option(t), 0, L_Option(t)/Endowment(t), R_Option(t)/Endowment(t)];
                    DG_SP = [DG_SP; update];
                end
            end
        end
        
     end
end
  
   

%% Run 2

if r2 > 0
    for r = 1
        
           % sub-101_task-ultimatum_run-0_raw.csv sub-102_task-ultimatum_run-1_raw.csv
        fname = fullfile(maindir,'logs',num2str(subj),sprintf('sub-%04d_task-ultimatum_run-%d_practice_raw.csv',subj,r)); % Psychopy taken out from Logs to make work for now.
        if exist(fname,'file')
            fid = fopen(fname,'r');
        else
            fprintf('sub-%d -- Let''s Make a Deal Game, Run %d: No data found.\n', subj, r+1)
            continue;
        end
        C = textscan(fid,repmat('%f',1,23),'Delimiter',',','HeaderLines',1,'EmptyValue', NaN);
        fclose(fid);
        
        
        % "Feedback" is the offer value (out of $20)
        
        Trial = C{1};
        Block = C{3};
        Endowment = C{4};
        response = C{17};
        response = round(response);
        L_Option = C{7};
        R_Option = C{8};
        
    end
    
    
   % Distribute games per responses
    
    % Block 1 = UG Prop
    % Block 2 = DG prop
    % Block 3 = UG Resp
    
     for t = 1:length(Endowment)
        
        if Block(t) == 4 % If DG Proposer
            if response(t) == 2 % If Left button selected
                if round(L_Option(t) > R_Option(t)) % Is Left button more?
                    update = [Trial(t), Endowment(t), L_Option(t), 1, L_Option(t)/Endowment(t), R_Option(t)/Endowment(t)];  % Trial number, Endowment, Choice, More, M_Option_Proportion, L_Option_Proportion
                    DG_FC = [DG_FC; update];
                else % Therefore it is less.
                    update = [Trial(t), Endowment(t), L_Option(t), 0, R_Option(t)/Endowment(t), L_Option(t)/Endowment(t)];
                    DG_FC = [DG_FC; update];
                end
            end
        end
        
        if Block(t) == 4
            if response(t) == 3 % Right option
                if round(R_Option(t) > L_Option(t)) % Is right option more?
                    update = [Trial(t), Endowment(t), R_Option(t), 1, R_Option(t)/Endowment(t), L_Option(t)/Endowment(t)];  % Trial number, Endowment, Choice, More, M_Option_Proportion, L_Option_Proportion
                    DG_FC = [DG_FC; update];
                else
                    update = [Trial(t), Endowment(t), R_Option(t), 0, L_Option(t)/Endowment(t), R_Option(t)/Endowment(t)];
                    DG_FC = [DG_FC; update];
                end
            end
        end
        
        
         if Block(t) == 3 % If DG Proposer
            if response(t) == 2 % If Left button selected
                if round(L_Option(t) > R_Option(t)) % Is Left button more?
                    update = [Trial(t), Endowment(t), L_Option(t), 1, L_Option(t)/Endowment(t), R_Option(t)/Endowment(t)];  % Trial number, Endowment, Choice, More, M_Option_Proportion, L_Option_Proportion
                    DG_SC = [DG_SC; update];
                else % Therefore it is less.
                    update = [Trial(t), Endowment(t), L_Option(t), 0, R_Option(t)/Endowment(t), L_Option(t)/Endowment(t)];
                    DG_SC = [DG_SC; update];
                end
            end
        end
        
        if Block(t) == 3
            if response(t) == 3 % Right option
                if round(R_Option(t) > L_Option(t)) % Is right option more?
                    update = [Trial(t), Endowment(t), R_Option(t), 1, R_Option(t)/Endowment(t), L_Option(t)/Endowment(t)];  % Trial number, Endowment, Choice, More, M_Option_Proportion, L_Option_Proportion
                    DG_SC = [DG_SC; update];
                else
                    update = [Trial(t), Endowment(t), R_Option(t), 0, L_Option(t)/Endowment(t), R_Option(t)/Endowment(t)];
                    DG_SC = [DG_SC; update];
                end
            end
        end
        
         if Block(t) == 2 % If DG Proposer
            if response(t) == 2 % If Left button selected
                if round(L_Option(t) > R_Option(t)) % Is Left button more?
                    update = [Trial(t), Endowment(t), L_Option(t), 1, L_Option(t)/Endowment(t), R_Option(t)/Endowment(t)];  % Trial number, Endowment, Choice, More, M_Option_Proportion, L_Option_Proportion
                    DG_FP = [DG_FP; update];
                else % Therefore it is less.
                    update = [Trial(t), Endowment(t), L_Option(t), 0, R_Option(t)/Endowment(t), L_Option(t)/Endowment(t)];
                    DG_FP = [DG_FP; update];
                end
            end
        end
        
        if Block(t) == 2
            if response(t) == 3 % Right option
                if round(R_Option(t) > L_Option(t)) % Is right option more?
                    update = [Trial(t), Endowment(t), R_Option(t), 1, R_Option(t)/Endowment(t), L_Option(t)/Endowment(t)];  % Trial number, Endowment, Choice, More, M_Option_Proportion, L_Option_Proportion
                    DG_FP = [DG_FP; update];
                else
                    update = [Trial(t), Endowment(t), R_Option(t), 0, L_Option(t)/Endowment(t), R_Option(t)/Endowment(t)];
                    DG_FP = [DG_FP; update];
                end
            end
        end
        
         if Block(t) == 1 % If DG Proposer
            if response(t) == 2 % If Left button selected
                if round(L_Option(t) > R_Option(t)) % Is Left button more?
                    update = [Trial(t), Endowment(t), L_Option(t), 1, L_Option(t)/Endowment(t), R_Option(t)/Endowment(t)];  % Trial number, Endowment, Choice, More, M_Option_Proportion, L_Option_Proportion
                    DG_SP = [DG_SP; update];
                else % Therefore it is less.
                    update = [Trial(t), Endowment(t), L_Option(t), 0, R_Option(t)/Endowment(t), L_Option(t)/Endowment(t)];
                    DG_SP = [DG_SP; update];
                end
            end
        end
        
        if Block(t) == 1
            if response(t) == 3 % Right option
                if round(R_Option(t) > L_Option(t)) % Is right option more?
                    update = [Trial(t), Endowment(t), R_Option(t), 1, R_Option(t)/Endowment(t), L_Option(t)/Endowment(t)];  % Trial number, Endowment, Choice, More, M_Option_Proportion, L_Option_Proportion
                    DG_SP = [DG_SP; update];
                else
                    update = [Trial(t), Endowment(t), R_Option(t), 0, L_Option(t)/Endowment(t), R_Option(t)/Endowment(t)];
                    DG_SP = [DG_SP; update];
                end
            end
        end
        
     end
end
  

%% DG Earnings

% DG_P is easy... simply subtract endowment from offer.

for ii = 1
    
    try
        
        a = size(DG_P);
        b = size(DG_P_2);
        a = a(1);
        b = b(1);
          
 
    if a>0 || b==0
        
        DG_P_Earnings = sum(DG_P(:,2) - DG_P(:,3));
    end
    
    if b>0 || a==0
        
        DG_P_Earnings = sum(DG_P_2(:,2) - DG_P_2(:,3));
        
    end
    
    if a>0 || b>0
        
        DG_P_Earnings = (sum(DG_P(:,2) - DG_P(:,3)) + sum(DG_P_2(:,2) - DG_P_2(:,3)))/2;
         
    end
    end
    
end

%% Save Proposer
  
concat = [DG_P;DG_P_2];
DG_P_Raw = array2table(concat(1:end,:),'VariableNames', {'Trial','Endowment','Choice','Decision','More_Prop','Less_Prop'});
output = ['Subject_' num2str(subj) '_DGP.csv'];
name = [output_folder,output];
writetable(DG_P_Raw, name); % Save as csv file

end

end
