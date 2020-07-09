classdef KF_tracker_OCC < handle
    % KF_tracker_OCC Kalman Filter track with occlusion model
    % A class to handle the data and processing of each track.
    % 
    % KF_tracker_OCC Properties:
    % dim - Sensor dimension
    % X - Position (x, y) and Velocity (vx, vy)
    % R - Measuremente Covariance
    % P - State Covariance
    % f - Model
    % event_rate - Events rate in us/event
    % prev_event - Timestamp of previous event
    % dt - Time difference between events
    % goodness - Represents the proximity of an event to the tracker location. High proximity gets a low score.
    % locked - Flag for locked state
    % active - Flag for active state
    % matched - Flag for matched state
    % dist_factor - Distance factor, assign events that are within a number of standard deviations from the location
    % occlusion - Occlusion flag
    % color - Color in RGB for drawing the boundary box
    
    properties
        dim
        X
        R
        P
        f = [1,0,1,0; 0,1,0,1; 0,0,1,0; 0,0,0,1];
        event_rate = 2e3;
        prev_event = 0;
        dt
        goodness = 1e6;
        locked = 0;
        active = 0;
        matched = 0;
        dist_factor = 8;
        occlusion = 0;
        color;
    end
    
    properties (Constant)
        H = [1,0,0,0; 0,1,0,0];                                             % dp/dx; p = pixel location of tracker
        border_thresh = 10;                                                 % deactivate trackers within # pixels of the border
        
        active_thresh = [0.2e5^2, 1e5^2];
        Q_free = diag([2, 2, 0.5, 0.5]);                                    % growth in state covariance with each assigned event. Large at first (since location and velocity are not locked)
        Q_active = diag([1, 1, 0.5, 0.5]);                                  % growth in state covariance with each assigned event. Smaller afer active
        Q_locked = diag([0.01, 0.01, 0.01, 0.01]);                          % growth in state covariance with each assigned event. Small afer locking
        
        alpha_free = 0.5;
        alpha_active = 0.05;
        alpha_locked = 0.5;
        
        min_size = 3;                                                       % minimum size of a tracker
        min_size_active = 3;
        
        size_factor = 2;
        max_vP = 150;
    end
    
    methods
        function obj = KF_tracker_OCC(dim, prev_event, color)
            obj.dim = dim;
            obj.f = [1,0,1,0; 0,1,0,1; 0,0,1,0; 0,0,0,1];
            obj.X = [dim(1)/2; dim(2)/2; 0; 0];                             % initial position (pixels) and velocity (pixels per event)
            obj.R = diag([100,100]);                                        % initial measurement covariance (pixels^2)
            obj.P = diag([25,25,25,25]);                                    % initial state covariance (pixels^2)
            obj.event_rate = 2e3;
            obj.locked = 0;
            obj.active = 0;
            obj.goodness = 1e6;
            obj.prev_event = prev_event;
            obj.color = color;
            obj.dist_factor = 3;
            obj.occlusion = 0;
        end
        
        function calculate_match(obj, xy)
            y = xy-obj.X(1:2);
            if obj.active == 0
                y = y/4;
            end
            
            obj.goodness = abs(y(1))/sqrt(obj.R(1,1)) ...
                + abs(y(2))/sqrt(obj.R(2,2));
            obj.matched = (abs(y(1)) < sqrt(obj.R(1,1)) * obj.dist_factor) ...
                && (abs(y(2)) < sqrt(obj.R(2,2)) * obj.dist_factor);
        end
        
        function update(obj,xy,ts)
            % update    Determine wich parameters should be used based on the state of the tracker
            obj.dt = (ts - obj.prev_event)/1e6;
            if obj.active == 1
                if obj.locked == 1
                    alpha = obj.alpha_locked;
                    Q = obj.Q_locked;
                else
                    alpha = obj.alpha_active;
                    Q = obj.Q_active;
                end
                alpha = exp(obj.dt*log(alpha));                             % make updates proportional to time
                Q = Q.*obj.dt;                                              % Q makes P grow proportional to time
            else
                alpha = obj.alpha_free;
                Q = obj.Q_free;
            end
            
            obj.f(1,3) = obj.dt;
            obj.f(2,4) = obj.dt;
            
            K = obj.updateX(xy, Q);                                         % update location of tracker
            if obj.occlusion == 0
                obj.updateR(xy, alpha);                                     % update size of tracker if no occlusion
            end
            
            obj.updateP(K);                                                 % update uncertainty
            obj.check_active(ts);
        end
        
        function K = updateX(obj, xy, Q)
            saved_state = obj.X;
            obj.X = obj.f*obj.X;
            y = xy - obj.X(1:2);                                            % error
            
            obj.P = obj.f*obj.P*obj.f' + Q;
            S = obj.H*obj.P*obj.H' + obj.R;
            K = obj.P*obj.H'/S;
            obj.X = obj.X + K*y;
            if obj.occlusion ~= 0                                           % keep previous velocities if occlusion occuring
                obj.X = [obj.X(1,1); obj.X(2,1); ...
                    saved_state(3,1); saved_state(4,1)];
            end
        end
        
        function updateR(obj, xy, alpha)
            y_res = xy - obj.X(1:2);                                        % error after updating X
            obj.R = alpha*obj.R + (1-alpha) * ...
                (y_res*y_res' + obj.H*obj.P*obj.H');
            
            if obj.active == 0
                obj.R(1,1) = max(obj.R(1,1), obj.min_size^2);
                obj.R(2,2) = max(obj.R(2,2), obj.min_size^2);
            else
                obj.R(1,1) = max(obj.R(1,1), obj.min_size_active^2);
                obj.R(2,2) = max(obj.R(2,2), obj.min_size_active^2);
            end
        end
        
        function updateP(obj,K)
            % finish Kalman update
            obj.P = (eye(4) - K*obj.H)*obj.P;
            obj.P(3,3) = min(obj.P(3,3), obj.max_vP);
            obj.P(4,4) = min(obj.P(4,4), obj.max_vP);
        end
        
        function check_active(obj, ts)
            obj.event_rate = 0.98*obj.event_rate ...
                + 0.02*(ts-obj.prev_event);
            obj.prev_event = ts;
            
            if (abs(obj.X(1)-obj.dim(1)) < obj.border_thresh) ...
                    || (abs(obj.X(1)) < obj.border_thresh) ...
                    || (abs(obj.X(2)-obj.dim(2)) < obj.border_thresh) ...
                    || (abs(obj.X(2)) < obj.border_thresh)
                obj.active = 0;
                obj.locked = 0;
                obj.X = [obj.dim(1)/2; obj.dim(2)/2; 0; 0];                 % initial position (pixels) and velocity (pixels per event)
                obj.R = diag([100,100]);                                    % initial measurement covariance (pixels^2)
                obj.P = diag([25,25,25,25]);                                % initial state covariance (pixels^2)
                obj.event_rate = 2e3;                                       % us/event
                obj.occlusion = 0;
            else
                if obj.active == 1
                    if obj.event_rate^2*(obj.R(1,1)*obj.R(2,2)) ...
                            < obj.active_thresh(2)
                        if ((obj.P(1,1) + obj.P(2,2)) < 4)
                            obj.locked = 1;
                        end
                    else
                        obj.active = 0;
                        obj.locked = 0;
                        obj.occlusion = 0;
                        obj.X = [obj.dim(1)/2; obj.dim(2)/2; 0; 0];
                        obj.R = diag([100,100]);
                        obj.P = diag([25,25,25,25]);
                        obj.event_rate = 2e3;
                    end
                else
                    if obj.event_rate^2*(obj.R(1,1)*obj.R(2,2)) ...
                            < obj.active_thresh(1)
                        obj.active = 1;
                        if ((obj.P(1,1) + obj.P(2,2)) < 4)
                            obj.locked = 1;
                        end
                    else
                        obj.active = 0;
                    end
                end
            end
        end
    end
end

