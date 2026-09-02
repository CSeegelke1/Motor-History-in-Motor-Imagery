% Title: Distinct contributions of motor imagery and execution to history-dependent biases in reaching (Experiment 2, A054)
% Authors: Seegelke, Heed 
% Script: Christian Seegelke 01/09/2026
% INPUT: A054_Data_complete_S02.mat
% OUTPUT: FIGURE 1B (Probe trajectories following execution and no movement trials of an exemplary participant (S02)
%========================================================================================================================

%%
clear; close all; clc

%% READ IN DATA
allSubjNr = [2]; % exemplary participant
data_AllSubj      = [];

datapath          = uigetdir(path, 'Select data directory');
% Define path for saving figures
figure_path = uigetdir(path, 'Select figure directory');

% load complete data
for i = allSubjNr
    if i <= 9
        load([datapath  '\' 'A054_Data_complete_S0' num2str(i)])
    else
        load([datapath  '\' 'A054_Data_complete_S' num2str(i)])
    end
    
    
    % raw data
    data_AllSubj      = [data_AllSubj data];
    
end


data = data_AllSubj;

%% Determine outliers as in R Script
% RT < 100ms & MT < 1000ms = outliers (both prime and probe)
for s = 1:size(data,2)
     for i = 1:size(data,1)
        if data(i,s).Error_Any == 0 && data(i,s).blocks_thisRepN ~= 0 && ( isnan(data(i,s).num_peak_vel_Prime) || data(i,s).num_peak_vel_Prime <=2 ) && data(i,s).num_peak_vel_Probe <=2 %no error trials and no practice trials and only smooth trajectories
            if data(i,s).RT_Prime < 100 || data(i,s).RT_Probe < 100 || data(i,s).MT_Prime > 1000 || data(i,s).MT_Probe > 1000
                data(i,s).Outlier = 1;
            else
                data(i,s).Outlier = 0;
            end
        else
            data(i,s).Outlier = NaN;
        end
     end
end
   

disp(['Outlier ', num2str(sum([data.Outlier],'omitna'))]) 


%% PROBE MOVEMENT: GET TARGET AND OBSTALCE POSITIONS FOR PLOTTING
% col 1: x-pos; col 2: y-pos; col 3: angle; col 4: radius major axix; col 5: radius minor axis

for i = 1:length(data)
    COORD_PROBE.start(1,1) = data(i).x_pos_start_location;
    COORD_PROBE.start(1,2) = data(i).y_pos_start_location;
    COORD_PROBE.start(1,3) = NaN; % no angle
    COORD_PROBE.start(1,4) = data(i).r1_start_location;
    COORD_PROBE.start(1,5) = data(i).r2_start_location; 
        if data(i).target_location_probe == 1
            COORD_PROBE.target(1,1) = data(i).x_pos_target_probe;
            COORD_PROBE.target(1,2) = data(i).y_pos_target_probe;
            COORD_PROBE.target(1,3) = data(i).deg_targets_probe;
            COORD_PROBE.target(1,4) = data(i).r1_target_location;
            COORD_PROBE.target(1,5) = data(i).r2_target_location;
        elseif data(i).target_location_probe == 2
            COORD_PROBE.target(2,1) = data(i).x_pos_target_probe;
            COORD_PROBE.target(2,2) = data(i).y_pos_target_probe;
            COORD_PROBE.target(2,3) = data(i).deg_targets_probe;
            COORD_PROBE.target(2,4) = data(i).r1_target_location;
            COORD_PROBE.target(2,5) = data(i).r2_target_location;
        elseif data(i).target_location_probe == 3
            COORD_PROBE.target(3,1) = data(i).x_pos_target_probe;
            COORD_PROBE.target(3,2) = data(i).y_pos_target_probe;
            COORD_PROBE.target(3,3) = data(i).deg_targets_probe;
            COORD_PROBE.target(3,4) = data(i).r1_target_location;
            COORD_PROBE.target(3,5) = data(i).r2_target_location;
        elseif data(i).target_location_probe == 4
            COORD_PROBE.target(4,1) = data(i).x_pos_target_probe;
            COORD_PROBE.target(4,2) = data(i).y_pos_target_probe;
            COORD_PROBE.target(4,3) = data(i).deg_targets_probe;
            COORD_PROBE.target(4,4) = data(i).r1_target_location;
            COORD_PROBE.target(4,5) = data(i).r2_target_location;
        elseif data(i).target_location_probe == 5
            COORD_PROBE.target(5,1) = data(i).x_pos_target_probe;
            COORD_PROBE.target(5,2) = data(i).y_pos_target_probe;
            COORD_PROBE.target(5,3) = data(i).deg_targets_probe;
            COORD_PROBE.target(5,4) = data(i).r1_target_location;
            COORD_PROBE.target(5,5) = data(i).r2_target_location;
        elseif data(i).target_location_probe == 6
            COORD_PROBE.target(6,1) = data(i).x_pos_target_probe;
            COORD_PROBE.target(6,2) = data(i).y_pos_target_probe;
            COORD_PROBE.target(6,3) = data(i).deg_targets_probe;
            COORD_PROBE.target(6,4) = data(i).r1_target_location;
            COORD_PROBE.target(6,5) = data(i).r2_target_location;
        elseif data(i).target_location_probe == 7
            COORD_PROBE.target(7,1) = data(i).x_pos_target_probe;
            COORD_PROBE.target(7,2) = data(i).y_pos_target_probe;
            COORD_PROBE.target(7,3) = data(i).deg_targets_probe;
            COORD_PROBE.target(7,4) = data(i).r1_target_location;
            COORD_PROBE.target(7,5) = data(i).r2_target_location;
        elseif data(i).target_location_probe == 8
            COORD_PROBE.target(8,1) = data(i).x_pos_target_probe;
            COORD_PROBE.target(8,2) = data(i).y_pos_target_probe;
            COORD_PROBE.target(8,3) = data(i).deg_targets_probe;
            COORD_PROBE.target(8,4) = data(i).r1_target_location;
            COORD_PROBE.target(8,5) = data(i).r2_target_location;
        else
        end
end

%obstacle coordinates
for i = 1:length(data)
        if data(i).target_location_probe == 1
            COORD_PROBE.obstacle(1,1) = data(i).x_pos_obs_location_probe;
            COORD_PROBE.obstacle(1,2) = data(i).y_pos_obs_location_probe;
            COORD_PROBE.obstacle(1,3) = data(i).deg_targets_probe;
            COORD_PROBE.obstacle(1,4) = data(i).r1_obstacle;
            COORD_PROBE.obstacle(1,5) = data(i).r2_obstacle;
        elseif data(i).target_location_probe == 2
            COORD_PROBE.obstacle(2,1) = data(i).x_pos_obs_location_probe;
            COORD_PROBE.obstacle(2,2) = data(i).y_pos_obs_location_probe;
            COORD_PROBE.obstacle(2,3) = data(i).deg_targets_probe;
            COORD_PROBE.obstacle(2,4) = data(i).r1_obstacle;
            COORD_PROBE.obstacle(2,5) = data(i).r2_obstacle;
        elseif data(i).target_location_probe == 3
            COORD_PROBE.obstacle(3,1) = data(i).x_pos_obs_location_probe;
            COORD_PROBE.obstacle(3,2) = data(i).y_pos_obs_location_probe;
            COORD_PROBE.obstacle(3,3) = data(i).deg_targets_probe;
            COORD_PROBE.obstacle(3,4) = data(i).r1_obstacle;
            COORD_PROBE.obstacle(3,5) = data(i).r2_obstacle;
        elseif data(i).target_location_probe == 4
            COORD_PROBE.obstacle(4,1) = data(i).x_pos_obs_location_probe;
            COORD_PROBE.obstacle(4,2) = data(i).y_pos_obs_location_probe;
            COORD_PROBE.obstacle(4,3) = data(i).deg_targets_probe;
            COORD_PROBE.obstacle(4,4) = data(i).r1_obstacle;
            COORD_PROBE.obstacle(4,5) = data(i).r2_obstacle;
        elseif data(i).target_location_probe == 5
            COORD_PROBE.obstacle(5,1) = data(i).x_pos_obs_location_probe;
            COORD_PROBE.obstacle(5,2) = data(i).y_pos_obs_location_probe;
            COORD_PROBE.obstacle(5,3) = data(i).deg_targets_probe;
            COORD_PROBE.obstacle(5,4) = data(i).r1_obstacle;
            COORD_PROBE.obstacle(5,5) = data(i).r2_obstacle;
        elseif data(i).target_location_probe == 6
            COORD_PROBE.obstacle(6,1) = data(i).x_pos_obs_location_probe;
            COORD_PROBE.obstacle(6,2) = data(i).y_pos_obs_location_probe;
            COORD_PROBE.obstacle(6,3) = data(i).deg_targets_probe;
            COORD_PROBE.obstacle(6,4) = data(i).r1_obstacle;
            COORD_PROBE.obstacle(6,5) = data(i).r2_obstacle;
        elseif data(i).target_location_probe == 7
            COORD_PROBE.obstacle(7,1) = data(i).x_pos_obs_location_probe;
            COORD_PROBE.obstacle(7,2) = data(i).y_pos_obs_location_probe;
            COORD_PROBE.obstacle(7,3) = data(i).deg_targets_probe;
            COORD_PROBE.obstacle(7,4) = data(i).r1_obstacle;
            COORD_PROBE.obstacle(7,5) = data(i).r2_obstacle;
        elseif data(i).target_location_probe == 8
            COORD_PROBE.obstacle(8,1) = data(i).x_pos_obs_location_probe;
            COORD_PROBE.obstacle(8,2) = data(i).y_pos_obs_location_probe;
            COORD_PROBE.obstacle(8,3) = data(i).deg_targets_probe;
            COORD_PROBE.obstacle(8,4) = data(i).r1_obstacle;
            COORD_PROBE.obstacle(8,5) = data(i).r2_obstacle;
        else
        end
end

%% PLOTTING PROBE MOVEMENTS: SINGLE TRIAL TRAJECTORIES FOR EACH SUBJECT OF PROBE MOVEMENTS
alpha = 0.95;
linewidth = 2;

for s = 1:size(data,2)
    
    Fig_Probe_Traj = figure('units','normalized','position',[0,0,0.54,0.96]); %h = suptitle(['S' num2str(data(1,s).subID) ' Probe movements']);h.Interpreter = 'none';
    h1 = subplot(2,2,1); title('with obstacle')
    p1 = get(h1,'position');
    set(gca,'linewidth',1,'FontSize',18)
    ellipse(COORD_PROBE.start(1,4),COORD_PROBE.start(1,5),0,COORD_PROBE.start(1,1),COORD_PROBE.start(1,2),'b'), hold on
    for i = 1:length(COORD_PROBE.target)
        ellipse(COORD_PROBE.target(i,4),COORD_PROBE.target(i,5),0,COORD_PROBE.target(i,1),COORD_PROBE.target(i,2),'k'), hold on
        ellipse(COORD_PROBE.obstacle(i,5),COORD_PROBE.obstacle(i,4),deg2rad(COORD_PROBE.obstacle(i,3)),COORD_PROBE.obstacle(i,1),COORD_PROBE.obstacle(i,2),'k'), hold on
    
        axis equal
    end

    for i = 1:size(data,1)
        if data(i,s).Error_Any == 0 && data(i,s).blocks_thisRepN ~= 0 && ( isnan(data(i,s).num_peak_vel_Prime) || data(i,s).num_peak_vel_Prime <=2 ) && data(i,s).num_peak_vel_Probe <=2 %no error trials and no practice trials and only smooth trajectories
            if strcmp(data(i,s).stop_signal_prime,'go') && strcmp(data(i,s).obstacle_probe,'yes')
                if strcmp(data(i,s).obstacle_prime,'no')
                    color = hex2rgb('#9C9C9C');
                elseif strcmp(data(i,s).obstacle_prime,'yes')
                    color = hex2rgb('#008000'); %dark green
                end
                obs_plot   = plot(data(i,s).mouse_probe_x_filt(data(i,s).mov_onset_Probe:data(i,s).mov_onset_Probe),data(i,s).mouse_probe_y_filt(data(i,s).mov_onset_Probe:data(i,s).mov_onset_Probe),'Linewidth',linewidth,'color', hex2rgb('#008000'));
                noobs_plot = plot(data(i,s).mouse_probe_x_filt(data(i,s).mov_onset_Probe:data(i,s).mov_onset_Probe),data(i,s).mouse_probe_y_filt(data(i,s).mov_onset_Probe:data(i,s).mov_onset_Probe),'Linewidth',linewidth,'color', hex2rgb('#9C9C9C'));
                patchline(data(i,s).mouse_probe_x_filt(data(i,s).mov_onset_Probe:data(i,s).frameNum_probe_time_cursor_in_target),data(i,s).mouse_probe_y_filt(data(i,s).mov_onset_Probe:data(i,s).frameNum_probe_time_cursor_in_target),'Linewidth',linewidth,'EdgeColor',color,'edgealpha',alpha), hold on
              
                xlim([-10 10])
                ylim([-10 10])
                
            end
        else
        end
    end

    lg_go_obs = legend([obs_plot noobs_plot], 'with obstacle','without obstacle');
    lg_go_obs.Position(1:2) = [.36 .88];
    title(lg_go_obs,'Prime')
    legend('boxoff')

    
    h2 = subplot(2,2,2);title('without obstacle')
    p2 = get(h2,'position');
    set(gca,'linewidth',1,'FontSize',18)
    ellipse(COORD_PROBE.start(1,4),COORD_PROBE.start(1,5),0,COORD_PROBE.start(1,1),COORD_PROBE.start(1,2),'b'), hold on
    for c = 1:length(COORD_PROBE.target)
        ellipse(COORD_PROBE.target(c,4),COORD_PROBE.target(c,5),0,COORD_PROBE.target(c,1),COORD_PROBE.target(c,2),'k'), hold on
        axis equal
    end
    
    for i = 1:size(data,1)
        if data(i,s).Error_Any == 0 && data(i,s).blocks_thisRepN ~= 0 && ( isnan(data(i,s).num_peak_vel_Prime) || data(i,s).num_peak_vel_Prime <=2 ) && data(i,s).num_peak_vel_Probe <=2 && data(i,s).Outlier ~= 1 %no error trials and no practice trials and only smooth trajectories, no outliers
            if strcmp(data(i,s).stop_signal_prime,'go') && strcmp(data(i,s).obstacle_probe,'no')
                if strcmp(data(i,s).obstacle_prime,'no')
                    color = hex2rgb('#9C9C9C'); %grey
                elseif strcmp(data(i,s).obstacle_prime,'yes')
                    color = hex2rgb('#008000'); %dark green
                end
                obs_plot   = plot(data(i,s).mouse_probe_x_filt(data(i,s).mov_onset_Probe:data(i,s).mov_onset_Probe),data(i,s).mouse_probe_y_filt(data(i,s).mov_onset_Probe:data(i,s).mov_onset_Probe),'Linewidth',linewidth,'color', hex2rgb('#008000'));
                noobs_plot = plot(data(i,s).mouse_probe_x_filt(data(i,s).mov_onset_Probe:data(i,s).mov_onset_Probe),data(i,s).mouse_probe_y_filt(data(i,s).mov_onset_Probe:data(i,s).mov_onset_Probe),'Linewidth',linewidth,'color', hex2rgb('#9C9C9C'));
                patchline(data(i,s).mouse_probe_x_filt(data(i,s).mov_onset_Probe:data(i,s).frameNum_probe_time_cursor_in_target),data(i,s).mouse_probe_y_filt(data(i,s).mov_onset_Probe:data(i,s).frameNum_probe_time_cursor_in_target),'Linewidth',linewidth,'EdgeColor',color,'edgealpha',alpha), hold on
        
                xlim([-10 10])
                ylim([-10 10])
                
            end
        else
        end
    end

    lg_go_noobs = legend([obs_plot noobs_plot], 'with obstacle','without obstacle');
    lg_go_noobs.Position(1:2) = [.80 .88];
    title(lg_go_noobs,'Prime')
    legend('boxoff')


    h3 = subplot(2,2,3); title('')
    p3 = get(h3,'position');
    set(gca,'linewidth',1,'FontSize',18)
    ellipse(COORD_PROBE.start(1,4),COORD_PROBE.start(1,5),0,COORD_PROBE.start(1,1),COORD_PROBE.start(1,2),'b'), hold on
    for i = 1:length(COORD_PROBE.target)
        ellipse(COORD_PROBE.target(i,4),COORD_PROBE.target(i,5),0,COORD_PROBE.target(i,1),COORD_PROBE.target(i,2),'k'), hold on
        ellipse(COORD_PROBE.obstacle(i,5),COORD_PROBE.obstacle(i,4),deg2rad(COORD_PROBE.obstacle(i,3)),COORD_PROBE.obstacle(i,1),COORD_PROBE.obstacle(i,2),'k'), hold on
    
        axis equal
    end

    for i = 1:size(data,1)
        if data(i,s).Error_Any == 0 && data(i,s).blocks_thisRepN ~= 0 && ( isnan(data(i,s).num_peak_vel_Prime) || data(i,s).num_peak_vel_Prime <=2 ) && data(i,s).num_peak_vel_Probe <=2 %no error trials and no practice trials and only smooth trajectories
            if strcmp(data(i,s).stop_signal_prime,'stop') && strcmp(data(i,s).obstacle_probe,'yes')
                if strcmp(data(i,s).obstacle_prime,'no')
                    color = hex2rgb('#9C9C9C');
                elseif strcmp(data(i,s).obstacle_prime,'yes')
                    color = hex2rgb('#808000'); %dark blue
                end
                patchline(data(i,s).mouse_probe_x_filt(data(i,s).mov_onset_Probe:data(i,s).frameNum_probe_time_cursor_in_target),data(i,s).mouse_probe_y_filt(data(i,s).mov_onset_Probe:data(i,s).frameNum_probe_time_cursor_in_target),'Linewidth',linewidth,'EdgeColor',color,'edgealpha',alpha), hold on
                obs_plot   = plot(data(i,s).mouse_probe_x_filt(data(i,s).mov_onset_Probe:data(i,s).mov_onset_Probe),data(i,s).mouse_probe_y_filt(data(i,s).mov_onset_Probe:data(i,s).mov_onset_Probe),'Linewidth',linewidth,'color', hex2rgb('#808000'));
                noobs_plot = plot(data(i,s).mouse_probe_x_filt(data(i,s).mov_onset_Probe:data(i,s).mov_onset_Probe),data(i,s).mouse_probe_y_filt(data(i,s).mov_onset_Probe:data(i,s).mov_onset_Probe),'Linewidth',linewidth,'color', hex2rgb('#9C9C9C'));
                xlim([-10 10])
                ylim([-10 10])
                
            end
        else
        end
    end
    
    lg_stop_obs = legend([obs_plot noobs_plot], 'with obstacle','without obstacle');
    lg_stop_obs.Position(1:2) = [.36 .4];
    title(lg_stop_obs,'Prime')
    legend('boxoff')



    h4 = subplot(2,2,4); title('')
    p4 = get(h4,'position');
    set(gca,'linewidth',1,'FontSize',18)
    ellipse(COORD_PROBE.start(1,4),COORD_PROBE.start(1,5),0,COORD_PROBE.start(1,1),COORD_PROBE.start(1,2),'b'), hold on
    for i = 1:length(COORD_PROBE.target)
        ellipse(COORD_PROBE.target(i,4),COORD_PROBE.target(i,5),0,COORD_PROBE.target(i,1),COORD_PROBE.target(i,2),'k'), hold on
        axis equal
    end

    for i = 1:size(data,1)
        if data(i,s).Error_Any == 0 && data(i,s).blocks_thisRepN ~= 0 && ( isnan(data(i,s).num_peak_vel_Prime) || data(i,s).num_peak_vel_Prime <=2 ) && data(i,s).num_peak_vel_Probe <=2 %no error trials and no practice trials and only smooth trajectories
            if strcmp(data(i,s).stop_signal_prime,'stop') && strcmp(data(i,s).obstacle_probe,'no')
                if strcmp(data(i,s).obstacle_prime,'no')
                    color = hex2rgb('#9C9C9C'); %grey
                elseif strcmp(data(i,s).obstacle_prime,'yes')
                    color = hex2rgb('#808000'); % dark blue
                end
                obs_plot   = plot(data(i,s).mouse_probe_x_filt(data(i,s).mov_onset_Probe:data(i,s).mov_onset_Probe),data(i,s).mouse_probe_y_filt(data(i,s).mov_onset_Probe:data(i,s).mov_onset_Probe),'Linewidth',linewidth,'color', hex2rgb('#808000'));
                noobs_plot = plot(data(i,s).mouse_probe_x_filt(data(i,s).mov_onset_Probe:data(i,s).mov_onset_Probe),data(i,s).mouse_probe_y_filt(data(i,s).mov_onset_Probe:data(i,s).mov_onset_Probe),'Linewidth',linewidth,'color', hex2rgb('#9C9C9C'));
                patchline(data(i,s).mouse_probe_x_filt(data(i,s).mov_onset_Probe:data(i,s).frameNum_probe_time_cursor_in_target),data(i,s).mouse_probe_y_filt(data(i,s).mov_onset_Probe:data(i,s).frameNum_probe_time_cursor_in_target),'Linewidth',linewidth,'EdgeColor',color,'edgealpha',alpha), hold on
            
                xlim([-10 10])
                ylim([-10 10])
                
            end
        else
        end
    end

    lg_stop_noobs = legend([obs_plot noobs_plot], 'with obstacle','without obstacle');
    lg_stop_noobs.Position(1:2) = [.80 .4];
    title(lg_stop_noobs,'Prime')
    legend('boxoff')

    
    % adjust X-label
    heightX=p2(2)+p2(4)-p4(2);
    h6=axes('position',[p3(1) p4(2) p4(3) heightX],'visible','off');
    v_label=xlabel('x-position [cm]','visible','on','FontSize',20);
    
    % adjust X-label
    heightX=p2(2)+p2(4)-p4(2);
    h6=axes('position',[p4(1) p4(2) p4(3) heightX],'visible','off');
    v_label=xlabel('x-position [cm]','visible','on','FontSize',20);
    
    % adjust Y-label
    heightY=p1(2)+p1(4)-p3(2);
    h5=axes('position',[p3(1) p3(2)-0.25 0.2 heightY],'visible','off');
    h_label=ylabel('y-position [cm]','visible','on','FontSize',20);
    h5_1=axes('position',[p3(1)-0.05 p3(2)-0.25 0.2 heightY],'visible','off');
    h_label1=ylabel('No Movement','visible','on','FontSize',25, 'FontWeight','bold');

    % adjust Y-label
    heightY=p1(2)+p1(4)-p3(2);
    h5=axes('position',[p3(1) p1(2)-0.25 0.2 heightY],'visible','off');
    h_label=ylabel('y-position [cm]','visible','on','FontSize',20);
    h6_1=axes('position',[p3(1)-0.05 p1(2)-0.25 0.2 heightY],'visible','off');
    h_label2=ylabel('Execution','visible','on','FontSize',25, 'FontWeight','bold');

    % export figure
    
    orient(Fig_Probe_Traj,'landscape')
    if data(1,s).subID <= 9
        %print(Fig_Probe_Traj,[figure_path, '\A054_ProbeTraj_S0', num2str(data(1,s).subID)],'-dpdf','-r0','-fillpage')
        exportgraphics(Fig_Probe_Traj,[figure_path, '\A054_ProbeTraj_S0', num2str(data(1,s).subID),'.pdf'],'ContentType','vector')
        saveas(Fig_Probe_Traj,[figure_path, '\A054_ProbeTraj_S0' num2str(data(1,s).subID) '.png'])
    else
        %print(Fig_Probe_Traj,[figure_path, '\A054_ProbeTraj_S', num2str(data(1,s).subID)],'-dpdf','-r0','-fillpage')
        exportgraphics(Fig_Probe_Traj,[figure_path, '\A054_ProbeTraj_S', num2str(data(1,s).subID),'.pdf'],'ContentType','vector')
        saveas(Fig_Probe_Traj,[figure_path, '\A054_ProbeTraj_S' num2str(data(1,s).subID) '.png'])
    end
end
