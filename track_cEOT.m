function [] = track_cEOT(filenamebin, verbose, show_frame)

close all;
%% load and prepare the data
fprintf('Loading bin file...\n\n');
TD = read_linux([filenamebin, '.bin']);

% fixes timestamp gaps
a = diff(TD.ts);
q = find(a>=2^16);
for qq = 1:length(q)
    TD.ts(q(qq):end) = TD.ts(q(qq):end) - 2^16;
end

TD = ImplementRefraction(TD, 10e3);                                         % 10ms refraction
TD = FilterTD(TD, 100e3);

%% init trackers
num_trackers = 16;
cc = jet(num_trackers);
for jj = 1:num_trackers
    tracker(jj) = KF_tracker_OCC([240,180], TD.ts(1), cc(jj,:));
end

%% init constants and variables
if show_frame == 1
    clf
    hold on
end

img_back = ones(max(TD.y),max(TD.x))*0.5;
frame_length = 1e6/24;
track_frame_length = 66e3;
track_frame_time = TD.ts(1);
frame_time = TD.ts(1);
last_index = 1;
timestep = 0.03;                                                            % in terms of events for occlusion detection

track_frame_num = 1;
k=1;
i=1;
v_count = 0;
vv_count = 0;

blank.x = zeros(1,num_trackers);
blank.y = zeros(1,num_trackers);
blank.w = zeros(1,num_trackers);
blank.h = zeros(1,num_trackers);
blank.locked = zeros(1,num_trackers);
blank.active = zeros(1,num_trackers);
blank.id = zeros(1,num_trackers);
blank.event_rate = zeros(1,num_trackers);
blank.eventRateArea = zeros(1,num_trackers);

fprintf('\nProcessing event by event:\n\n');
%% event by event processing
for ii = 1:length(TD.ts)
    
    %% show sanity check output
    if verbose == 1
        fprintf('Event at timestamp: %d\n', TD.ts(ii));
    else
%         if v_count == 50
%             fprintf('\n.');
%             vv_count = vv_count + 1;
%             if vv_count == v_count
%                 vv_count = 0;
%             end
%             v_count = 0;
%         else
%             if v_count == vv_count
%                 fprintf('*');
%             else
%                 fprintf('.');
%             end
%         end
%         v_count = v_count + 1;
        pattern_ev = "..................................................";
%         fprintf("%d - %d\n", TD.y(ii), mod(TD.y(ii), 50));
        pattern_ev{1}(mod(TD.y(ii), 50)+1) = '*';
        pattern_ev = char(pattern_ev);
        fprintf("%s\n", pattern_ev);
    end
    
    %% compute distance for each tracker and check whether it is within range
    for jj = 1:num_trackers
        tracker(jj).calculate_match([TD.x(ii);TD.y(ii)]);
    end
    
    %% find which tracker to assign the event to
    if sum([tracker(logical([tracker.locked])).matched]>0)                  % if at least one locked tracker is matched
        min_goodness = inf;
        for jj = 1:num_trackers
            if(tracker(jj).locked ...
                    && tracker(jj).matched) ...
                    && (tracker(jj).goodness < min_goodness)
                
                min_goodness = tracker(jj).goodness;
                update_tracker = jj;
            end
        end
    else
        if sum([tracker(logical([tracker.active])).matched]>0)              % if at least one active tracker is matched
            min_goodness = inf;
            for jj = 1:num_trackers
                if(tracker(jj).active ...
                        && tracker(jj).matched) ...
                        && (tracker(jj).goodness < min_goodness)
                    
                    min_goodness = tracker(jj).goodness;
                    update_tracker = jj;
                end
            end
        else
            if sum([tracker.matched]>0)
                [~, update_tracker] = min([tracker.goodness]);
            else
                update_tracker = 0;
            end
        end
    end
    
    %% check if multiple locked trackers match event
    check_occlusion = zeros(1,16);
    if sum([tracker(logical([tracker.locked])).matched]>0) > 1              % if at least two locked tracker are matched
        for jj = 1:num_trackers
            check_occlusion(jj) = tracker(jj).locked ...
                && tracker(jj).matched;
        end
    end
    
    %% check for overlaps for up to 2 time steps among locked trackers that match to same event
    overlap_trackers = [];
    if any(check_occlusion)                                                 % check for overlaps
        matched_trackers = find(check_occlusion);
        for jj = 1:length(matched_trackers)
            for kk = (jj+1):length(matched_trackers)
                if abs(tracker(matched_trackers(jj)).X(3)-tracker(matched_trackers(kk)).X(3))>50 ...
                        || (abs(tracker(matched_trackers(jj)).X(3)-tracker(matched_trackers(kk)).X(3))>20 && ((tracker(matched_trackers(jj)).X(3)>0 && tracker(matched_trackers(kk)).X(3)<0) ...
                        || (tracker(matched_trackers(jj)).X(3)<0 && tracker(matched_trackers(kk)).X(3)>0)))
                    
                    if sum(diag(tracker(matched_trackers(jj)).P))<1.5 ...
                            && sum(diag(tracker(matched_trackers(kk)).P))<1.5
                        
                        %x:left, y:bottom
                        p1_x = tracker(matched_trackers(jj)).X(1,1)-sqrt(tracker(matched_trackers(jj)).R(1,1))*tracker(matched_trackers(jj)).size_factor;
                        p1_y = tracker(matched_trackers(jj)).X(2,1)+sqrt(tracker(matched_trackers(jj)).R(2,2))*tracker(matched_trackers(jj)).size_factor;
                        p1_xvelocity = tracker(matched_trackers(jj)).X(3,1);
                        p1_yvelocity = tracker(matched_trackers(jj)).X(4,1);
                        p1_w = sqrt(tracker(matched_trackers(jj)).R(1,1))*2*tracker(matched_trackers(jj)).size_factor;
                        p1_h = sqrt(tracker(matched_trackers(jj)).R(2,2))*2*tracker(matched_trackers(jj)).size_factor;
                        p2_x = tracker(matched_trackers(kk)).X(1,1)-sqrt(tracker(matched_trackers(kk)).R(1,1))*tracker(matched_trackers(kk)).size_factor;
                        p2_y = tracker(matched_trackers(kk)).X(2,1)+sqrt(tracker(matched_trackers(kk)).R(2,2))*tracker(matched_trackers(kk)).size_factor;
                        p2_xvelocity = tracker(matched_trackers(kk)).X(3,1);
                        p2_yvelocity = tracker(matched_trackers(kk)).X(4,1);
                        p2_w = sqrt(tracker(matched_trackers(kk)).R(1,1))*2*tracker(matched_trackers(kk)).size_factor;
                        p2_h = sqrt(tracker(matched_trackers(kk)).R(2,2))*2*tracker(matched_trackers(kk)).size_factor;
                        
                        p1_step1 = [p1_x+p1_xvelocity*timestep p1_y+p1_yvelocity*timestep p1_w p1_h];
                        p1_step2 = [p1_x+p1_xvelocity*timestep*2 p1_y+p1_yvelocity*timestep*2 p1_w p1_h];
                        p2_step1 = [p2_x+p2_xvelocity*timestep p2_y+p2_yvelocity*timestep p2_w p2_h];
                        p2_step2 = [p2_x+p2_xvelocity*timestep*2 p2_y+p2_yvelocity*timestep*2 p2_w p2_h];
                        
                        if (rectdist(p1_step1, p2_step1)==0 ...
                                || rectdist(p1_step2, p2_step2)==0)
                            
                            tracker(matched_trackers(jj)).occlusion = matched_trackers(kk);
                            tracker(matched_trackers(kk)).occlusion = matched_trackers(jj);
                            if ~ismember(matched_trackers(jj), overlap_trackers)
                                overlap_trackers = cat(2, overlap_trackers, matched_trackers(jj));
                            end
                            if ~ismember(matched_trackers(kk), overlap_trackers)
                                overlap_trackers = cat(2, overlap_trackers, matched_trackers(kk));
                            end
                        end
                    end
                end
            end
        end
    end
    
    %% implement the update
    if any(overlap_trackers)                                                % if overlapping, update tracker based on previous velocities
        for aa = 1:length(overlap_trackers)
            tracker(overlap_trackers(aa)).update([TD.x(ii);TD.y(ii)], TD.ts(ii));
        end
    else
        if update_tracker>0                                                 % no multiple matching trackers
            if (tracker(update_tracker).occlusion~=0 ...
                    && tracker(update_tracker).locked==1)                   % check if previously occluded
                
                if tracker(tracker(update_tracker).occlusion).locked==1
                    tracker1 = update_tracker;
                    tracker2 = tracker(update_tracker).occlusion;
                    
                    p1_x = tracker(tracker1).X(1,1)-sqrt(tracker(tracker1).R(1,1))*tracker(tracker1).size_factor;
                    p1_y = tracker(tracker1).X(2,1)+sqrt(tracker(tracker1).R(2,2))*tracker(tracker1).size_factor;
                    p1_xvelocity = tracker(tracker1).X(3,1);
                    p1_yvelocity = tracker(tracker1).X(4,1);
                    p1_w = sqrt(tracker(tracker1).R(1,1))*2*tracker(tracker1).size_factor;
                    p1_h = sqrt(tracker(tracker1).R(2,2))*2*tracker(tracker1).size_factor;
                    p2_x = tracker(tracker2).X(1,1)-sqrt(tracker(tracker2).R(1,1))*tracker(tracker2).size_factor;
                    p2_y = tracker(tracker2).X(2,1)+sqrt(tracker(tracker2).R(2,2))*tracker(tracker2).size_factor;
                    p2_xvelocity = tracker(tracker2).X(3,1);
                    p2_yvelocity = tracker(tracker2).X(4,1);
                    p2_w = sqrt(tracker(tracker2).R(1,1))*2*tracker(tracker2).size_factor;
                    p2_h = sqrt(tracker(tracker2).R(2,2))*2*tracker(tracker2).size_factor;
                    
                    p1_step0 = [p1_x p1_y p1_w p1_h];
                    p1_step1 = [p1_x+p1_xvelocity*timestep p1_y+p1_yvelocity*timestep p1_w p1_h];
                    p1_step2 = [p1_x+p1_xvelocity*timestep*2 p1_y+p1_yvelocity*timestep*2 p1_w p1_h];
                    p2_step0 = [p2_x p2_y p2_w p2_h];
                    p2_step1 = [p2_x+p2_xvelocity*timestep p2_y+p2_yvelocity*timestep p2_w p2_h];
                    p2_step2 = [p2_x+p2_xvelocity*timestep*2 p2_y+p2_yvelocity*timestep*2 p2_w p2_h];
                    if rectdist(p1_step0, p2_step0)>0 ...
                            && rectdist(p1_step1, p2_step1)>0 ...
                            && rectdist(p1_step2, p2_step2)>0               % if no longer occluding
                        
                        tracker(tracker1).occlusion = 0;
                        tracker(tracker2).occlusion = 0;
                    end
                else
                    tracker(tracker1).occlusion = 0;
                    tracker(tracker2).occlusion = 0;
                end
            end
            tracker(update_tracker).update([TD.x(ii);TD.y(ii)], TD.ts(ii)); % update normally
        end
    end
    
    %% display the data as a frame
    if (TD.ts(ii)>frame_time)
        current = blank;
        %% periodic cleanup
        for jj = 1:num_trackers
            tracker(jj).update(tracker(jj).X(1:2), TD.ts(ii));
            if tracker(jj).active == 0
                tracker(jj).R = tracker(jj).R + diag([20,20]);
            end
            
            if tracker(jj).locked == 1
                for kk = (jj+1):num_trackers
                    if tracker(kk).locked == 1
                        % compute the vertices of the overlapping area of the rectangles
                        left    = max(tracker(jj).X(1) - sqrt(tracker(jj).R(1,1))*tracker(jj).size_factor, tracker(kk).X(1) - sqrt(tracker(kk).R(1,1))*tracker(kk).size_factor);
                        right   = min(tracker(jj).X(1) + sqrt(tracker(jj).R(1,1))*tracker(jj).size_factor, tracker(kk).X(1) + sqrt(tracker(kk).R(1,1))*tracker(kk).size_factor);
                        top     = max(tracker(jj).X(2) - sqrt(tracker(jj).R(2,2))*tracker(jj).size_factor, tracker(kk).X(2) - sqrt(tracker(kk).R(2,2))*tracker(kk).size_factor);
                        bottom  = min(tracker(jj).X(2) + sqrt(tracker(jj).R(2,2))*tracker(jj).size_factor, tracker(kk).X(2) + sqrt(tracker(kk).R(2,2))*tracker(kk).size_factor);
                        
                        % compute overlapping area (intersection)
                        if (left<=(right-5) && top <= (bottom-5)) && (norm(tracker(jj).X(3:4) - tracker(kk).X(3:4))<20)
                            if ((tracker(jj).X(1) > left) && (tracker(jj).X(1) < right)) ...
                                    || ((tracker(kk).X(1) > left) && (tracker(kk).X(1) < right))
                                
                                tracker(jj).R(1,1) =  mean([sqrt(tracker(jj).R(1,1)), sqrt(tracker(kk).R(1,1))])^2;
                            else
                                tracker(jj).R(1,1) =  (sqrt(tracker(jj).R(1,1)) + sqrt(tracker(kk).R(1,1)))^2;
                            end
                            
                            if ((tracker(jj).X(2) > top) && (tracker(jj).X(2) < bottom)) ...
                                    || ((tracker(kk).X(2) > top) && (tracker(kk).X(2) < bottom))
                                
                                tracker(jj).R(2,2) =  mean([sqrt(tracker(jj).R(2,2)), sqrt(tracker(kk).R(2,2))])^2;
                            else
                                tracker(jj).R(2,2) =  (sqrt(tracker(jj).R(2,2)) + sqrt(tracker(kk).R(2,2)))^2;
                            end
                            
                            tracker(jj).X = (tracker(jj).X+tracker(kk).X)/2;
                            
                            tracker(kk) =  KF_tracker_OCC([240,180], TD.ts(1), tracker(kk).color);
                        else                                                % if tracker boxes are near and moving in same direction but not overlapping
                            if abs(right-left)<=5 && norm(tracker(jj).X(3:4) - tracker(kk).X(3:4))<20 && ((tracker(jj).X(3)<0 && tracker(kk).X(3)<0) ...
                                    || (tracker(jj).X(3)>0 && tracker(kk).X(3)>0))
                                
                                tracker(jj).R(1,1) =  (sqrt(tracker(jj).R(1,1)) + sqrt(tracker(kk).R(1,1)))^2;
                                tracker(jj).R(2,2) =  mean([sqrt(tracker(jj).R(2,2)), sqrt(tracker(kk).R(2,2))])^2;
                                tracker(jj).X =  (tracker(jj).X+tracker(kk).X)/2;
                                tracker(kk) =  KF_tracker_OCC([240,180], TD.ts(1), tracker(kk).color);
                            end
                        end
                    end
                end
            end
        end
        
        %% show frame
        frame_time = frame_time + frame_length;
        
        indices = sub2ind(size(img_back), TD.y(last_index:ii), TD.x(last_index:ii));
        img = img_back;
        img(indices) = 2-TD.p(last_index:ii);
        last_index = ii;
        if show_frame == 1
            clf
            imshow(img);
            for jj = 1:num_trackers
                if tracker(jj).active == 1
                    rectangle('Position', ...
                        [tracker(jj).X(1)-sqrt(tracker(jj).R(1,1))*tracker(jj).size_factor, ...
                        tracker(jj).X(2)-sqrt(tracker(jj).R(2,2))*tracker(jj).size_factor, ...
                        sqrt(tracker(jj).R(1,1))*2*tracker(jj).size_factor, ...
                        sqrt(tracker(jj).R(2,2))*2*tracker(jj).size_factor], ...
                        'EdgeColor', tracker(jj).color);
                    if tracker(jj).locked == 1
                        rectangle('Position', ...
                            [tracker(jj).X(1)-sqrt(tracker(jj).R(1,1))*tracker(jj).size_factor, ...
                            tracker(jj).X(2)-sqrt(tracker(jj).R(2,2))*tracker(jj).size_factor, ...
                            5, 5], ...
                            'FaceColor', 'g');
                    end
                end
            end
            drawnow
        end
        
        if show_frame == 1
            vid(k) = getframe;
        end
        k = k+1;
    end
    
    %% For output file
    if (TD.ts(ii) > track_frame_time)
        if(TD.ts(ii) > (track_frame_time+track_frame_length))
            while(TD.ts(ii) > (track_frame_time+track_frame_length))
                tracks{track_frame_num} = blank;
                track_frame_num = track_frame_num + 1;
                track_frame_time = track_frame_time + track_frame_length;
            end
        end
        for jj = 1:num_trackers
            current.x(jj) = tracker(jj).X(1);
            current.y(jj) = tracker(jj).X(2);
            current.w(jj) = sqrt(tracker(jj).R(1,1))*2*tracker(jj).size_factor;
            current.h(jj) = sqrt(tracker(jj).R(2,2))*2*tracker(jj).size_factor;
            current.locked(jj) = tracker(jj).locked;
            current.active(jj) = tracker(jj).active;
            current.event_rate(jj) = tracker(jj).event_rate;
            current.eventRateArea(jj) = tracker(jj).event_rate^2*(tracker(jj).R(1,1)*tracker(jj).R(2,2));
            if (track_frame_num ~= 1 ...
                    && ~isequal(tracks{track_frame_num-1}, blank))
                
                if (TD.ts(ii) < (track_frame_time+track_frame_length) ...
                        && tracker(jj).locked == 1 ...
                        && tracks{track_frame_num-1}.locked(jj) == 1)
                    
                    current.id(jj) = tracks{track_frame_num-1}.id(jj);
                else
                    if (tracker(jj).locked == 1)
                        current.id(jj) = i;
                        i = i+1;
                    end
                end
            else
                if (tracker(jj).locked == 1)
                    current.id(jj) = i;
                    i = i+1;
                end
            end
        end
        tracks{track_frame_num} = current;
        track_frame_num = track_frame_num + 1;
        track_frame_time = track_frame_time + track_frame_length;
    end
end

if show_frame == 1
    SaveMovie(vid, [filenamebin '_Vid_Occ'])
end

fprintf('\n\nGenerating output file...\n\n');
%% Generate output file
fileID = fopen([filenamebin, '_full.txt'], 'w');
ini_comments = ['#This text file contains annotation data for recordings in: Recordings\\', filenamebin,'.bin\n', ...
    '#The corresponding picture of the recording site is at: Image\\', filenamebin,'.png\n', ...
    '#The annotation file is stored at: NUS Tracker Annotations\\', filenamebin,'.txt\n', ...
    '#Comments: Annotated with NUS tracker code\n', ...
    '#The recordings are not annotated by: Anyone\n', ...
    '#Sensor Dimensions- Height = 180 Pixels Width = 240 Pixels\n', ...
    '#LEGEND: 1-Car, 2-Bus, 3-Pedestrian, 4-Bike, 5-Truck, 6-Unknown\n', ...
    '#Time(us),x-Location,y-Location,x-size,y-size,track-num,class\n'];
fprintf(fileID, ini_comments);

for i = 1:length(tracks)
    if ~isequal(tracks{i}, blank)
        for index = 1:num_trackers
            if tracks{i}.locked(index) == 1
                fprintf(fileID, '%.0f, %d, %d, %u, %u, %u, 6,\n', TD.ts(1)+((i-1)*track_frame_length), ...
                    round( tracks{i}.x(index) ), ...
                    round( tracks{i}.y(index) ), ...
                    round( tracks{i}.w(index) ), ...
                    round( tracks{i}.h(index) ), ...
                    tracks{i}.id(index));
            end
        end
    end
end
fclose(fileID);