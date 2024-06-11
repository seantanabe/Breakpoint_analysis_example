
clear all 
cd D:\20240610_Montreal_anes\breakpoint_example
load('Breakpoint_onset_end_example_dt_simu.mat') % load attatched data

% preparing simulated dt which resembles 
% EMG gamma flexor activity during propofol induction
psd_tmp = zeros(length(t_dt), 1)
psd_tmp(1:100) = 0.3; 
psd_tmp(100:150) = 0.3:((1-0.3)/50):1; 
psd_tmp(150:350) = 1; 
psd_tmp(350:400) = 1:(-1/50):0; 
rng(0)
psd_tmp = psd_tmp + rand(length(t_dt),1)*0.3; 
psd_tmp = imgaussfilt(psd_tmp,20);
psd_tmp = rescale(psd_tmp); 
% figure; plot(psd_tmp)

% brk_tmp: 
% time index of breakpoints determined by structural breakpoin analysis

% c_i is 7, for gamma band extensor/flexor muscle  
% c_i is the signal type (EEG/EMG, PSD band)
% Since certain signal types increase/decrease during induction

% psd_tmp: power spectral density

% t_dt: time (s)

min_diff = 0.08;

psd_brk_diff = diff(psd_tmp(brk_tmp));

% find onset/end

if any(c_i == [1 2 3]) % Cz delta, theta, alpha
    for fold_i = 1 % find increase onset/end
        tmp1   = find(psd_brk_diff > min_diff);
        tmp2_a = diff(tmp1');
        tmp2_b = find([tmp2_a inf]>1);
        tmp2_c = diff([0 tmp2_b]); %length of the sequences
        tmp2_d = cumsum(tmp2_c) + 1; % endpoints of the sequences
        tmp2_e = tmp2_d - tmp2_c ; %onset

        tmp_plot      = [tmp2_e; tmp2_d]';
        tmp_plot(:,2) = tmp_plot(:,2) - 1;
        tmp_plot1 = NaN(size(tmp_plot));
        tmp_plot1(:,1)= tmp1(tmp_plot(:,1));
        tmp_plot1(:,2)= tmp1(tmp_plot(:,2));
        tmp_plot1(:,2)= tmp_plot1(:,2) + 1;

        [~,tmp_indbrk] = max(psd_tmp(brk_tmp(tmp_plot1(:,2))) - psd_tmp(brk_tmp(tmp_plot1(:,1))));
        tmp_plot_brk =  [brk_tmp(tmp_plot1(tmp_indbrk,1)) brk_tmp(tmp_plot1(tmp_indbrk,2))];
        tmp_plot_brk_pos = tmp_plot_brk;
    end
elseif any(c_i == [4 6 7]) % Cz beta, flex/ext gamma
    for fold_i = 1 % find increase onset/end
        tmp1   = find(psd_brk_diff > min_diff);
        tmp2_a = diff(tmp1');
        tmp2_b = find([tmp2_a inf]>1);
        tmp2_c = diff([0 tmp2_b]); %length of the sequences
        tmp2_d = cumsum(tmp2_c) + 1; % endpoints of the sequences
        tmp2_e = tmp2_d - tmp2_c ; %onset

        tmp_plot      = [tmp2_e; tmp2_d]';
        tmp_plot(:,2) = tmp_plot(:,2) - 1;
        tmp_plot1 = NaN(size(tmp_plot));
        tmp_plot1(:,1)= tmp1(tmp_plot(:,1));
        tmp_plot1(:,2)= tmp1(tmp_plot(:,2));
        tmp_plot1(:,2)= tmp_plot1(:,2) + 1;

        [~,tmp_indbrk] = max(psd_tmp(brk_tmp(tmp_plot1(:,2))) - psd_tmp(brk_tmp(tmp_plot1(:,1))));
        tmp_plot_brk =  [brk_tmp(tmp_plot1(tmp_indbrk,1)) brk_tmp(tmp_plot1(tmp_indbrk,2))];
        tmp_plot_brk_pos = tmp_plot_brk;
    end
    for fold_i = 1 % find decrease onset/end
        tmp1   = find(psd_brk_diff < -min_diff);
        tmp2_a = diff(tmp1');
        tmp2_b = find([tmp2_a inf]>1);
        tmp2_c = diff([0 tmp2_b]); %length of the sequences
        tmp2_d = cumsum(tmp2_c) + 1; % endpoints of the sequences
        tmp2_e = tmp2_d - tmp2_c ; %onset

        tmp_plot      = [tmp2_e; tmp2_d]';
        tmp_plot(:,2) = tmp_plot(:,2) - 1;
        tmp_plot1 = NaN(size(tmp_plot));
        tmp_plot1(:,1)= tmp1(tmp_plot(:,1));
        tmp_plot1(:,2)= tmp1(tmp_plot(:,2));
        tmp_plot1(:,2)= tmp_plot1(:,2) + 1;

        [~,tmp_indbrk] = min(psd_tmp(brk_tmp(tmp_plot1(:,2))) - psd_tmp(brk_tmp(tmp_plot1(:,1))));
        tmp_plot_brk =  [brk_tmp(tmp_plot1(tmp_indbrk,1)) brk_tmp(tmp_plot1(tmp_indbrk,2))];
        tmp_plot_brk_neg = tmp_plot_brk;
    end
elseif any(c_i == [5]) % Cz gamma
    for fold_i = 1 % find decrease onset/end
        tmp1   = find(psd_brk_diff < -min_diff);
        tmp2_a = diff(tmp1');
        tmp2_b = find([tmp2_a inf]>1);
        tmp2_c = diff([0 tmp2_b]); %length of the sequences
        tmp2_d = cumsum(tmp2_c) + 1; % endpoints of the sequences
        tmp2_e = tmp2_d - tmp2_c ; %onset

        tmp_plot      = [tmp2_e; tmp2_d]';
        tmp_plot(:,2) = tmp_plot(:,2) - 1;
        tmp_plot1 = NaN(size(tmp_plot));
        tmp_plot1(:,1)= tmp1(tmp_plot(:,1));
        tmp_plot1(:,2)= tmp1(tmp_plot(:,2));
        tmp_plot1(:,2)= tmp_plot1(:,2) + 1;

        [~,tmp_indbrk] = min(psd_tmp(brk_tmp(tmp_plot1(:,2))) - psd_tmp(brk_tmp(tmp_plot1(:,1))));
        tmp_plot_brk =  [brk_tmp(tmp_plot1(tmp_indbrk,1)) brk_tmp(tmp_plot1(tmp_indbrk,2))];
        tmp_plot_brk_neg = tmp_plot_brk;
    end
end

figure; plot(t_dt, psd_tmp,'k') % ext gamma
hold on ; scatter(t_dt(tmp_plot_brk_pos(1)),psd_tmp(tmp_plot_brk_pos(1)),800,'+r')
hold on ; scatter(t_dt(tmp_plot_brk_pos(2)),psd_tmp(tmp_plot_brk_pos(2)),800,'Xr')
hold on ; scatter(t_dt(tmp_plot_brk_neg(1)),psd_tmp(tmp_plot_brk_neg(1)),800,'+b')
hold on ; scatter(t_dt(tmp_plot_brk_neg(2)),psd_tmp(tmp_plot_brk_neg(2)),800,'Xb'); hold off
ylabel('PSD'); xlabel('time (s)')
box off
set(gcf,'color','w'); set(gca,'TickDir','out');



