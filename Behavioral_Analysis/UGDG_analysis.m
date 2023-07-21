clear all
close all
clc

% This code imports subject earnings (across conditions) and their
% individual difference measures and analyzes the data. 

% Daniel Sazhin
% 02/15/2023
% DVS Lab
% Temple University

inputdir = fullfile(pwd,'Behavioral_Results\');
%% subjects

values = [10103, 10348, 10374, 10391, 10402, 10436, 10460, 10478, 10496, 10512, 10531, 10541, 10570, 10581, 10608, 12001, 12002, 12003, 12004];


%% Raw DG

subjects = values;
DG_P_Earnings = [];
DG_P_Offers=[];
DG_P_Prop = [];

for ii = 1:length(subjects)
    
    name = ['Subject_' num2str(subjects(ii)) '_DGP.csv'];
    input = [inputdir,name];
    O = readtable(input);
    DG_Part = [];
    for kk = 1:length(O.Trial)
        if O.Decision(kk) == 1
            save = O.More_Prop(kk);
        else
            save = O.Less_Prop(kk);
        end
        DG_Part = [DG_Part; save];
    end
    DG_P_Prop = [DG_P_Prop; mean(DG_Part)];
    DG_P_Offers = [DG_P_Offers; mean(O.Choice)]; % AMOUNT OFFERED!!!
    DG_P_Earnings = [DG_P_Earnings; sum(O.Endowment)-sum(O.Choice)]; % AMOUNT SAVED for self.
    
end


%% Raw results UG-P

Final_save_2 = [];
UG_P_2 = [];

Final_save = [];
UG_P = [];
%UG_P_Total = [];
Subjects = [];
Subjects_2 = [];
UG_P_Offers = [];
UG_P_Offers_2 = [];
UG_P_Earnings  = [];
UG_P_Earnings_2 = [];

subjects = values;
Final_Subjects =[];
UG_P_Prop = [];

for jj = 1:length(subjects)
    save_value = [];
    
    name = ['Subject_' num2str(subjects(jj)) '_UGP.csv'];
    input = [inputdir,name];
    T = readtable(input);
    UG_P_Offers = [UG_P_Offers; mean(T.Choice)]; % AMOUNT OFFERED!!!
    UG_P_Offers_Prop = [UG_P_Offers; mean(T.Choice)];
    
    UG_Part = [];
    
    for kk = 1:length(T.Trial)
        if T.Decision(kk) == 1
            save = T.More_Prop(kk);
        else
            save = T.Less_Prop(kk);
        end
        UG_Part = [UG_Part; save];
        
    end
    UG_P_Prop = [UG_P_Prop; mean(UG_Part)];
    UG_P_Earnings = [UG_P_Earnings; sum(T.Endowment)-sum(T.Choice)]; % AMOUNT SAVED for self, uncorrected for rejections.
  
end

UG_P_Prop_2 = [];

for jj = 1:length(subjects)
    save_value = [];
    name = ['Subject_' num2str(subjects(jj)) '_UGP2.csv'];
    input = [inputdir,name];
    S = readtable(input);
    
    UG_Part = [];
    
    for kk = 1:length(T.Trial)
        if T.Decision(kk) == 1
            save = T.More_Prop(kk);
        else
            save = T.Less_Prop(kk);
        end
        UG_Part = [UG_Part; save];
        
    end
    
    UG_P_Prop_2 = [UG_P_Prop_2; mean(UG_Part)];
    UG_P_Offers_2 = [UG_P_Offers_2;mean(S.Choice)]; % AMOUNT OFFERED!!!
    UG_P_Earnings_2 = [UG_P_Earnings_2; sum(S.Endowment)-sum(S.Choice)]; % AMOUNT SAVED for self, uncorrected for rejections.
    
end

UG_P_Raw_Props = ((UG_P_Prop+UG_P_Prop_2)/2);
UG_P_Raw_Offers = (UG_P_Offers+UG_P_Offers_2/2); % Sum offers
UG_P_Raw_Earnings = UG_P_Earnings + UG_P_Earnings_2; % Uncorrected
%% Plot proportions chosen

data = [mean(DG_P_Prop), mean(UG_P_Raw_Props)]; 
x = linspace(1,2,2);
fig = figure;
x1 = bar(x(1),data(1));
x1.LineWidth= 2.5;
hold on
x2 = bar(x(2),data(2));
x2.LineWidth= 2.5;
colormap('jet')

ax = gca;
ax.FontSize = 12;
box off
xlabel ('Tasks', 'FontSize', 16);
ylabel  ('Offers', 'FontSize', 16);
title ('Mean Endowment Offered By Participants')
set(gcf,'color','w');
set(gca, 'XTick', 1:2, 'XTickLabels', {'DG-P', 'UG-P'})

hold on

% Standard Error

B1Er = std(DG_P_Prop) / sqrt(length(DG_P_Prop));
B2Er = std(UG_P_Raw_Props) / sqrt(length(UG_P_Raw_Props));


err = [B1Er,B2Er] * 2;

er = errorbar(x,data,err); 
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 2.5;
hold off

saveas(gcf,'Bar_Props.tif')

%% Combine tasks in model.

% Import results.
% Combine into single matrix
% Add +1 for DG and -1 for UG
% Generate interaction effects. 

% Import data:

missing_data=[];
save_X = [];
save_Y = [];
Cor_Endow_DeltaP=[];
Cor_Endow_Intercept=[];
Cor_DeltaP_Intercept= [];
Main_Effects = [];
subject_regressor = [];

for ii = 1:length(values)
    S = [];
    T = [];

        try
            name = ['Subject_' num2str(values(ii)) '_UGP.csv'];
            input = [inputdir,name];
            S = readtable(input);
            name = ['Subject_' num2str(values(ii)) '_UGP2.csv'];
            input = [inputdir,name];
            T = readtable(input);
            subject_UGP=[S;T];
            
            name = ['Subject_' num2str(values(ii)) '_DGP.csv'];
            input = [inputdir,name];
            U = readtable(input);
            subject_DGP = [U];
            
        catch
            missing_data = [missing_data, values(ii)];
        end
        
        % Generate regressors
        
        UG_X = [subject_UGP.Endowment,(subject_UGP.More_Prop-subject_UGP.Less_Prop)]; % Generate endowment and delta P.
        DG_X = [subject_DGP.Endowment,(subject_DGP.More_Prop-subject_DGP.Less_Prop)];
        X = zscore([UG_X;DG_X]); % z scores demeans
        
        % Add Task regressor and Interaction of Task/Endowment, Task/Delta
        % P. DG is 1 and UG is -1.
        
        [N,M]= size(UG_X);
        UG_mat = 0*(ones(N,1)); % UG is coded as 0
        
        [N,M]= size(DG_X);
        DG_mat = ones(N,1); % DG is coded as 0
        
        task_regressor = [UG_mat; DG_mat];
        
        [N,M]= size(task_regressor);
        
        regessor_mat = ones(N,1);
        regessor_mat = regessor_mat*values(ii);

        subject_regressor = [subject_regressor; regessor_mat];
        
        % Generate interaction of task regressor and endowment and delta p.
        
        endow_int = task_regressor.*X(:,1);
        deltaP_int = task_regressor.*X(:,2);
        %mean_int = task_regressor.*X(:,3);
        
        % Concatenate into regressor matrix.
        Main_Effects_part = [task_regressor, X];
        X = [task_regressor, X, endow_int, deltaP_int]; 
        Y = [subject_UGP.Decision;subject_DGP.Decision];
        
        % Loop over all of the participants.
        save_X = [save_X; X];
        save_Y = [save_Y; Y];
        Main_Effects = [Main_Effects; Main_Effects_part];
        
%         % Save correlations of regressors across participants.

%         [B,DEV,STATS] = glmfit(X,Y,'binomial','link','logit');
%         coeff=STATS.coeffcorr;
%         Cor_Endow_DeltaP=[Cor_Endow_DeltaP; coeff(2,3)];
%         Cor_Endow_Intercept=[Cor_Endow_Intercept; coeff(1,2)];
%         Cor_DeltaP_Intercept= [Cor_DeltaP_Intercept; coeff(1,3)];
end


% Run logistic regression.

[B,DEV,STATS] = glmfit(save_X,save_Y,'binomial','link','logit');
%Y=categorical(save_Y)
%[B,dev,stats] = mnrfit(save_X,Y);
se=STATS.se; 
data = [B(1), B(2), B(3), B(4), B(5), B(6)]; % Intercept, Task, Endowment, Delta P, Endow Interaction, Delta P interaction.
x = linspace(1,6,6);
figure
bar(x,data)
ax = gca;
ax.FontSize = 12;
box off
xlabel ('Regressors', 'FontSize', 16);
ylabel  ('Z-standardized Beta Weight', 'FontSize', 16);
set(gcf,'color','w');
set(gca, 'XTick', 1:6, 'XTickLabels', {'Intercept', 'Task', 'Endowment', 'Delta P', 'Task*Endow', 'Task*DeltaP'})
title('Does Endowment and Offers Predict Proposer More/Less Choices?')

hold on

% Standard Error

B1Er = se(1);
B2Er = se(2);
B3Er = se(3);
B4Er = se(4);
B5Er = se(5);
B6Er = se(6);
%B7Er = se(7);
%B8Er = se(8);


err = [B1Er,B2Er,B3Er,B4Er,B5Er,B6Er] * 2;

er = errorbar(x,data,err); 
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 1;
hold off

saveas(gcf,'Bar_proposer.png')

%% Mixed effects model


tb1 = array2table(Main_Effects(1:end,:),'VariableNames', {'x1','x2','x3'}); %x1 is Task, x2 is Endowment, and x3 is Delta P. 
tb2 = array2table(save_Y,'VariableNames', {'y'});
tb3 = array2table(subject_regressor,'VariableNames', {'subject'});
mixed_effects = [tb2, tb3, tb1];

% y ~ x1 + x2 + x3 + intA + intB + (x1 | subject) + (x2 | subject) + (x3 | subject)

% Where y is "predictor", x1 is Task, x2 is Endowment, and x3 is Delta P. 


random_intercept_model = fitlme(mixed_effects,'y ~ x1 + x2 + x3 + x1:x2 + x1:x3 + (1 | subject)');
null_intercept = fitlme(mixed_effects,'y ~ 1 + (1 | subject)');

[TABLE,SIMINFO] = compare(null_intercept, random_intercept_model)

random_slope_model = fitlme(mixed_effects,'y ~ x1 + x2 + x3 + x1:x2 + x1:x3 + (x1 | subject) + (x2 | subject) + (x3 | subject)');
null_slope = fitlme(mixed_effects,'y ~ 1 + (x1 | subject) + (x2 | subject) + (x3 | subject)');

[TABLE,SIMINFO] = compare(null_slope, random_slope_model)

