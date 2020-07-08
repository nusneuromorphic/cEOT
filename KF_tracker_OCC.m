classdef KF_tracker_OCC < handle
    %KF_TRACKER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dim
        X
        R
        P
        event_rate = 2e3; %us/event
        prev_event = 0;
        locked = 0;
        active = 0;
        goodness = 1e6; %how good is the event match
        matched = 0;
        f = [1,0,1,0; 0,1,0,1; 0,0,1,0; 0,0,0,1]; %model
        color;
        dist_factor = 8; % assign events that are within 10 standard deviations of the location
        dt
        occlusion = 0;
    end
    
    properties (Constant)
        H = [1,0,0,0; 0,1,0,0]; %dp/dx ; p= pixel location of tracker
        border_thresh = 10; %deactivate trackers within 10 pixels of the border
        
        active_thresh = [0.2e5^2, 1e5^2];
        %active_thresh = [0.2e5^2, 0.9e5^2]; %testing active_thresh
        Q_free = diag([2, 2, 0.5, 0.5]); %growth in state covariance with each assigned event. Large at first (since location and velocity are not locked). Small afer locking.
        Q_active = diag([1, 1, 0.5, 0.5]); %growth in state covariance with each assigned event. Large at first (since location and velocity are not locked). Small afer locking.
        Q_locked = diag([0.01, 0.01, 0.01, 0.01]); %growth in state covariance with each assigned event. Large at first (since location and velocity are not locked). Small afer locking.
        
        
%         alpha_free = 0.8; %alpha can change from small (when determining size) to large when locked (to roughly maintain size)
%         alpha_active = 0.96; %alpha can change from small (when determining size) to large when locked (to roughly maintain size)
%         alpha_locked = 0.9999; %alpha can change from small (when determining size) to large when locked (to roughly maintain size)
        alpha_free = 0.5; %alpha can change from small (when determining size) to large when locked (to roughly maintain size)
        alpha_active = 0.05; %alpha can change from small (when determining size) to large when locked (to roughly maintain size)
        alpha_locked = 0.5; %alpha can change from small (when determining size) to large when locked (to roughly maintain size)
        
        %min_size = 3; %minimum size of a tracker
        min_size = 3; %minimum size of a tracker
        min_size_active = 3;
        %max_size = 20; %minimum size of a tracker
        size_factor = 2;
        max_vP = 150;
    end
    
    methods
        function obj = KF_tracker_OCC(dim, prev_event, color)
            obj.dim = dim;
            obj.f = [1,0,1,0; 0,1,0,1; 0,0,1,0; 0,0,0,1]; %model
            %obj.X = [rand*dim(1); rand*dim(2); 0; 0]; %intial position (pixels) and velocity (pixels per event)
            obj.X = [dim(1)/2; dim(2)/2; 0; 0];
            obj.R = diag([100,100]); %initial measurement covariance (pixels^2)
            obj.P = diag([25,25,25,25]); %initial state covariance (pixels^2)
            obj.event_rate = 2e3; %us/event
            %obj.event_rate = 1e6/24; %testing event_rate
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
            obj.goodness = abs(y(1))/sqrt(obj.R(1,1)) + abs(y(2))/sqrt(obj.R(2,2));
            obj.matched = (abs(y(1)) < sqrt(obj.R(1,1))*obj.dist_factor) && (abs(y(2)) < sqrt(obj.R(2,2))*obj.dist_factor);
        end
        
        function update(obj,xy,ts)
            %determine wich parameters should be used based on the state of
            %the tracker
            obj.dt = (ts - obj.prev_event)/1e6;
            if obj.active == 1
                if obj.locked == 1
                    alpha = obj.alpha_locked;
                    Q = obj.Q_locked;
                else
                    alpha = obj.alpha_active;
                    Q = obj.Q_active;
                end
                alpha = exp(obj.dt*log(alpha)); %make updates proportional to time
%                alpha = 1 - (1-alpha)*obj.dt; %make updates proportional to time
                Q = Q.*obj.dt; %Q makes P grow proportional to time
            else
                alpha = obj.alpha_free;
                Q = obj.Q_free;
            end

            obj.f(1,3) = obj.dt; %means X is pixels/sec
            obj.f(2,4) = obj.dt;
%             if obj.active == 1
%                 obj.f(1,3) = (ts - obj.prev_event)/1e6; %means X is pixels/sec
%                 obj.f(2,4) = (ts - obj.prev_event)/1e6;
%             else
%                 obj.f(1,3) = 0; %means X is pixels/sec
%                 obj.f(2,4) = 0;
%             end
            
            K = obj.updateX(xy, Q); %update location of tracker
            if obj.occlusion==0
                obj.updateR(xy, alpha); %update size of tracker if no occlusion
            end
            obj.updateP(K); %update uncertainty
            obj.check_active(ts);
        end
        
        function K = updateX(obj, xy, Q)
            saved_state = obj.X;
            obj.X = obj.f*obj.X;
            y = xy-obj.X(1:2); %error
            
            obj.P = obj.f*obj.P*obj.f' + Q;
            S = obj.H*obj.P*obj.H'+ obj.R;
            K = obj.P*obj.H'/S;
            obj.X = obj.X + K*y;
            if obj.occlusion~=0 %keep previous velocities if occlusion occuring
                obj.X = [obj.X(1,1); obj.X(2,1); saved_state(3,1); saved_state(4,1)];
            end
        end
        
        function updateR(obj, xy, alpha)
            y_res = xy-obj.X(1:2); %error after updating X
            obj.R = alpha*obj.R + (1-alpha)*(y_res*y_res' + obj.H*obj.P*obj.H');
            if obj.active == 0
                obj.R(1,1) = max(obj.R(1,1), obj.min_size^2);
                obj.R(2,2) = max(obj.R(2,2), obj.min_size^2);
            else
                obj.R(1,1) = max(obj.R(1,1), obj.min_size_active^2);
                obj.R(2,2) = max(obj.R(2,2), obj.min_size_active^2);
            end
%             else
%                 if min(obj.R(1,1), obj.R(2,2))  < obj.min_size_inactive^2
%                     obj.R = diag(obj.min_size_inactive^2, obj.min_size_inactive^2);
%                 end
                
%                 obj.R(1,1) = min(obj.R(1,1), obj.max_size^2);
%                 obj.R(2,2) = min(obj.R(2,2), obj.max_size^2);
            %end
            
%             obj.R(1,1) = min(obj.R(1,1), obj.max_size^2);
%             obj.R(2,2) = min(obj.R(2,2), obj.max_size^2);
            
%             if obj.active == 0
%                 obj.R = obj.R + 5*eye(2);
%             else
%                 if obj.locked == 0
%                     obj.R = obj.R + 0.01*eye(2);
%                 end
%             end
            
        end
        
        function updateP(obj,K)
            %finish Kalman update
            obj.P = (eye(4)-K*obj.H)*obj.P;
            obj.P(3,3) = min(obj.P(3,3), obj.max_vP);
            obj.P(4,4) = min(obj.P(4,4), obj.max_vP);
        end
        
        function check_active(obj, ts)
            obj.event_rate = 0.98*obj.event_rate + 0.02*(ts-obj.prev_event);
            obj.prev_event = ts;
            
            if (abs(obj.X(1)-obj.dim(1))< obj.border_thresh) || (abs(obj.X(1))< obj.border_thresh) || (abs(obj.X(2)-obj.dim(2))< obj.border_thresh) || (abs(obj.X(2))< obj.border_thresh)
                obj.active = 0;
                obj.locked = 0;
                %obj.X = [rand*obj.dim(1); rand*obj.dim(2); 0; 0]; %intial position (pixels) and velocity (pixels per event)
                obj.X = [obj.dim(1)/2; obj.dim(2)/2; 0; 0];
                obj.R = diag([100,100]); %initial measurement covariance (pixels^2)
                obj.P = diag([25,25,25,25]); %initial state covariance (pixels^2)
                obj.event_rate = 2e3; %us/event
                %obj.event_rate = 1e6/24; %testing event_rate
                obj.occlusion = 0;
            else
                if obj.active == 1
                    if obj.event_rate^2*(obj.R(1,1)*obj.R(2,2))<obj.active_thresh(2)
                        %obj.active = 1;
                        %if ((obj.P(1,1) + obj.P(2,2))<4) && ((obj.P(3,3) + obj.P(4,4))<200)
                        if ((obj.P(1,1) + obj.P(2,2))<4)
                            obj.locked = 1;
                        end
                    else
                        %fprintf('locked change\n')
                        obj.active = 0;
                        obj.locked = 0;
                        %obj.R(1,1) = 900;
                        %obj.R(2,2) = 900;
                        obj.occlusion = 0;
                        obj.X = [obj.dim(1)/2; obj.dim(2)/2; 0; 0];
                        obj.R = diag([100,100]); %initial measurement covariance (pixels^2)
                        obj.P = diag([25,25,25,25]); %initial state covariance (pixels^2)
                        obj.event_rate = 2e3; %us/event
                        %obj.event_rate = 1e6/24; %testing event_rate
                    end
                else
                    if obj.event_rate^2*(obj.R(1,1)*obj.R(2,2))<obj.active_thresh(1)
                        obj.active = 1;
                        %if ((obj.P(1,1) + obj.P(2,2))<4) && ((obj.P(3,3) + obj.P(4,4))<30)
                        if ((obj.P(1,1) + obj.P(2,2))<4)
                            obj.locked = 1;
                        end
                        %fprintf('locked change\n')
                    else
                        obj.active = 0;
                    end
                end
            end
        end
    end
    
end

