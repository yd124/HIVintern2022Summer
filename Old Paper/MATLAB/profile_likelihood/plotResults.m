cd("~/Documents/GitHub/HIVintern2022Summer/Old Paper/MATLAB/profile_likelihood")
% cd('~/HIVintern2022Summer/Old Paper/MATLAB/profile_likelihood')
data = readtable('../../data/data1/logvRNA.csv');
patients = unique(data.patient,'stable');

names = ["b0" "bi" "k" "dlt" "p" "d" "tau"];

for i=1:1
% for i=1:length(patients)
    patient = string(patients(i));

    cd("~/Documents/GitHub/HIVintern2022Summer/Old Paper/MATLAB/profile_likelihood/"+patient);
    LSEs = table2array(readtable('LSE22.csv'));
    sampled_params = table2array(readtable('sampled_params2.csv'));
    M = table2array(readtable(patient+'.csv'));
    
    num_params = size(LSEs,1);

    tiledlayout(4,2);
    for ip = 1:num_params
        nexttile
        P = sampled_params(ip,:);
        L = LSEs(ip,:);
        fig = scatter(P(L>0),L(L>0),'filled');
%         fig = plot(P(L>0),L(L>0));
        
        xlabel(names(ip));
        ylabel("LSE");
        hold on

        for j = 2:2 %size(M,1)
            fig = scatter(M(j,1+ip),M(j,1),'filled');
            hold on
        end
        legend('','x2','Location','northeast')

%         legend('','x1','x2','Location','northeast')

%         legend('','x1','x2','x3', 'x4','x5','x6','Location','northwest')
        saveas(fig,patient+'31.jpeg')
        hold off
    end 
end 


% LSEs = table2array(readtable('LSE1.csv'));
% sampled_params = table2array(readtable('sampled_params.csv'));
% % original_params = table2array(readtable('original_params.csv'));
% results = table2array(readtable('../new_work/results.csv'));
% original_params = table2array(readtable('original_params.csv'));
% 
% 
% xx_error = original_params(1);
% xx = original_params(2:end);
% 
% num_params = size(LSEs,1);
% 
% names = ["b0" "bi" "k" "dlt" "p" "d" "tau"];
% 
% tiledlayout(4,2);
% for ip = 1:num_params
%     nexttile
%     p1 = sampled_params(ip,:);
%     a1 = LSEs(ip,:);
% %     scatter(p1(a1>0),a1(a1>0),'filled');
%     plot(p1(a1>0),a1(a1>0));
%     xlabel(names(ip));
%     ylabel("LSE");
%     hold on
%     for i = 1:size(results,1)
%         scatter(results(i,1+ip),results(i,1),'filled');
%         hold on
%     end
%     legend('','x1','x2','x3', 'x4','x5','x6','Location','northwest')
% 
%     hold off
% end 