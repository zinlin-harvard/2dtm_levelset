classdef LevelSet
    %LevelSet: class consisting of a grid g, a level-set function phi, and 
    %          methods for updating and modifying the level-set function
    %          
    %    TODO: implement reinitialization, implement 3rd-order Runge-Kutta
    
    properties
        % type: grid2d
        g
        
        % array on grid.x, grid.y
        phi;
        phi_pad;
        x;
        y;
        
        % "int"
        reinit_step
        
        %
        max_num_steps = 100;
        
        % "double"
        alpha
        smooth_eps;
%         smooth_eps
%         dt
%         dtau
    end
    
    % methods
    %  LevelSet(varargin)
    %  update_fw_euler(obj, vn, t)
    %  dt_cfl
    %  update_rk3tvd
    %  grad_phi_hj_weno(obj, vn)
    
    methods
        % args: Nx, Ny, x0, y0, dx, dy, alpha, init_phi_fx
        %    or g, alpha, init_phi_fx
        %    or obj_in, phi (copy obj_in, with new phi)
        % init_phi_fx = function
        function obj = LevelSet(varargin) % can't overload constructor in Matlab...
            if (nargin == 3)
                obj.g = varargin{1};
                obj.alpha = varargin{2};
                init_phi_fx = varargin{3};
                obj.x = reshape(obj.g.x, obj.g.Nx, obj.g.Ny).';
                obj.y = reshape(obj.g.y, obj.g.Nx, obj.g.Ny).';
                obj.phi = init_phi_fx(obj.x, obj.y);
            elseif (nargin == 8)
                obj.g = grid2d(varargin{1}, varargin{2}, varargin{3}, ...
                    varargin{4}, varargin{5}, varargin{6});
                obj.alpha = varargin{7};
                init_phi_fx = varargin{8};
                obj.x = reshape(obj.g.x, obj.g.Nx, obj.g.Ny).';
                obj.y = reshape(obj.g.y, obj.g.Nx, obj.g.Ny).';
                obj.phi = init_phi_fx(obj.x, obj.y);
            elseif (nargin == 2)
                obj_in = varargin{1};
                obj.g = obj_in.g;
                obj.x = reshape(obj.g.x, obj.g.Nx, obj.g.Ny).';
                obj.y = reshape(obj.g.y, obj.g.Nx, obj.g.Ny).';
                obj.alpha = obj_in.alpha;
                obj.phi = varargin{2};
            end            
            obj.smooth_eps = 1.5 * (obj.g.dx + obj.g.dy) / 2;
            
            obj = obj.update_phi_pad;
        end
        
        % normal velocity vn, total time t
        % 1st-order time (forward difference), using 5th-order spatial accuracy
        function [obj, grad_phi] = update_fw_euler(obj, vx, vy, t)
            [dt, n_steps] = obj.dt_cfl(vx, vy, t);
            if(n_steps > obj.max_num_steps)
                error('number of time steps, %g, is greater than max number, %g', n_steps, obj.max_num_steps);
            end
            for i=1:n_steps
                [phi_x, phi_y] = grad_phi_hj_weno(obj, vx, vy);
                obj.phi = obj.phi - dt * (vx .* phi_x + vy .* phi_y);
            end
            obj = obj.update_phi_pad(); % need to update phi_pad when update phi
        end
        
        
        % CFL condition
        function [dt, num_steps] = dt_cfl(obj, vx, vy, t)
            dt = obj.alpha * min(obj.g.dx, obj.g.dy) / max(max(abs(vx(:))), max(abs(vy(:))));
            num_steps = ceil(t / dt);
            dt = t / num_steps;
        end
        
        % 5th-order spatial accuracy, 3rd-order time (Runge-Kutta)
        % v = vn or [vx vy]
        function obj = update_rk3tvd(obj, vx, vy, t)
            [dt, n_steps] = obj.dt_cfl(vx, vy, t);
            if(n_steps > obj.max_num_steps)
                error('number of time steps, %g, is greater than max number, %g', n_steps, obj.max_num_steps);
            end
            for i=1:n_steps
                [phi_x, phi_y] = grad_phi_hj_weno(obj, vx, vy);
                phi1 = obj.phi - dt * (vx .* phi_x + vy .* phi_y);
                obj1 = LevelSet(obj, phi1);
                [phi_x1, phi_y1] = grad_phi_hj_weno(obj1, vx, vy);
                phi2 = phi1 - dt * (vx .* phi_x1 + vy .* phi_y1);
                phi_half = 0.75 * obj.phi + 0.25 * phi2;
                obj_half = LevelSet(obj, phi_half);
                [phi_x_half, phi_y_half] = grad_phi_hj_weno(obj_half, vx, vy);
                phi_three_halves = phi_half - dt * (vx .* phi_x_half + vy .* phi_y_half);
                obj.phi = obj.phi / 3 + 2 * phi_three_halves / 3;
%                 obj.phi = obj.phi - dt * vn;
            end
            obj = obj.update_phi_pad(); % need to update phi_pad when update phi
        end
        
        function obj = reinit_rk3tvd(obj)
            for i=1:obj.reinit_step
                phi1 = reinit_kernal(obj);
                obj1 = LevelSet(obj, phi1);
                phi2 = reinit_kernal(obj1);
                phi_half = 0.75 * obj.phi + 0.25 * phi2;
                obj_half = LevelSet(obj, phi_half);
                phi_three_halves = reinit_kernal(obj_half);
                obj.phi = obj.phi / 3 + 2 * phi_three_halves / 3;
            end
            obj = obj.update_phi_pad(); % need to update phi_pad when update phi
        end
        
        % This is the kernal for reinitilization
        function [new_phi] = reinit_kernal(obj)
            dtau = 0.1 * min(obj.g.dx, obj.g.dy);
            epsilon = 2.0 * min(obj.g.dx, obj.g.dy);

            cent_ind = (4:size(obj.phi_pad,1)-3).';
            
            diff_x_pad = (obj.phi_pad(cent_ind,2:end) - obj.phi_pad(cent_ind,1:end-1)) / obj.g.dx;
            v1p = diff_x_pad(:,6:end);
            v2p = diff_x_pad(:,5:end-1);
            v3p = diff_x_pad(:,4:end-2);
            v4p = diff_x_pad(:,3:end-3);
            v5p = diff_x_pad(:,2:end-4);
            v5m = diff_x_pad(:,1:end-5);
            der_x_p = LevelSet.hj_weno_dphi(v1p, v2p, v3p, v4p, v5p);
            der_x_m = LevelSet.hj_weno_dphi(v5m, v5p, v4p, v3p, v2p);
            
            diff_y_pad = (obj.phi_pad(2:end,cent_ind) - obj.phi_pad(1:end-1,cent_ind)) / obj.g.dy;
            v1p = diff_y_pad(6:end,:);
            v2p = diff_y_pad(5:end-1,:);
            v3p = diff_y_pad(4:end-2,:);
            v4p = diff_y_pad(3:end-3,:);
            v5p = diff_y_pad(2:end-4,:);
            v5m = diff_y_pad(1:end-5,:);
            der_y_p = LevelSet.hj_weno_dphi(v1p, v2p, v3p, v4p, v5p);
            der_y_m = LevelSet.hj_weno_dphi(v5m, v5p, v4p, v3p, v2p);
            
            g_x_squared = max( max(der_x_m,0).^2, min(der_x_p,0).^2 ) .* (obj.phi>0) ...
                + max( min(der_x_m,0).^2, max(der_x_p,0).^2 ) .* (obj.phi<=0);
            g_y_squared = max( max(der_y_m,0).^2, min(der_y_p,0).^2 ) .* (obj.phi>0) ...
                + max( min(der_y_m,0).^2, max(der_y_p,0).^2 ) .* (obj.phi<=0);

            smoothing_factor = obj.phi ./ sqrt(obj.phi.^2 + epsilon.^2) ;
            new_phi = obj.phi - dtau * smoothing_factor .* (sqrt(g_x_squared + g_y_squared) - 1.0);
        end
        
        function [phi_x, phi_y] = grad_phi_hj_weno(obj, vx, vy)
            % v_plus = [v1p v2p v3p v4p v5p]
            % v_minus = [v2p v3p v4p v5p v5m]
            cent_ind = (4:size(obj.phi_pad,1)-3).';
            
            diff_x_pad = (obj.phi_pad(cent_ind,2:end) - obj.phi_pad(cent_ind,1:end-1)) / obj.g.dx;
            v1p = diff_x_pad(:,6:end);
            v2p = diff_x_pad(:,5:end-1);
            v3p = diff_x_pad(:,4:end-2);
            v4p = diff_x_pad(:,3:end-3);
            v5p = diff_x_pad(:,2:end-4);
            v5m = diff_x_pad(:,1:end-5);
            phi_x_p = LevelSet.hj_weno_dphi(v1p, v2p, v3p, v4p, v5p);
            %             phi_x_m = LevelSet.hj_weno_dphi(v2p, v3p, v4p, v5p, v5m);
            phi_x_m = LevelSet.hj_weno_dphi(v5m, v5p, v4p, v3p, v2p);
            
            diff_y_pad = (obj.phi_pad(2:end,cent_ind) - obj.phi_pad(1:end-1,cent_ind)) / obj.g.dy;
            v1p = diff_y_pad(6:end,:);
            v2p = diff_y_pad(5:end-1,:);
            v3p = diff_y_pad(4:end-2,:);
            v4p = diff_y_pad(3:end-3,:);
            v5p = diff_y_pad(2:end-4,:);
            v5m = diff_y_pad(1:end-5,:);
            phi_y_p = LevelSet.hj_weno_dphi(v1p, v2p, v3p, v4p, v5p);
            %             phi_y_m = LevelSet.hj_weno_dphi(v2p, v3p, v4p, v5p, v5m);
            phi_y_m = LevelSet.hj_weno_dphi(v5m, v5p, v4p, v3p, v2p);
            
            %             vn = v;
            %             nx_grad_phi = obj.phi_pad(cent_ind,5:end-2) - obj.phi_pad(cent_ind,3:end-4);
            %             ny_grad_phi = obj.phi_pad(5:end-2,cent_ind) - obj.phi_pad(3:end-4,cent_ind);
            %             vx = sign(vn) .* nx_grad_phi; % technically vx*grad_phi*mag
            %             vy = sign(vn) .* ny_grad_phi;
            
            eps = 1e-10; % only use phi_x is nx is nonzero, same for phi_y (Fan W.)
            phi_x = phi_x_p .* (vx<0) + phi_x_m .* (vx>=0);
            phi_y = phi_y_p .* (vy<0) + phi_y_m .* (vy>=0);
            
            %             grad_phi = sqrt(phi_x.*phi_x + phi_y.*phi_y);
        end

        function[val] = ls_interp(obj, x, y)
            F = griddedInterpolant(obj.x.', obj.y.', obj.phi.');
            val = F(x(:), y(:));
        end
        
        % "surface-integral" approximated by a volume integral with a
        % smoothed delta function (so f_vol should be defined on the grid
        % g). We use the smoothing function on pgs. 15,16 of Fedkiw & Osher
        function[val] = surface_integral(obj, fn_vol)
            val = sum(sum(fn_vol .* obj.delta_function_x));
        end
        
        % ignore scaling factors of 1/dx, 1/dy that are canceled in any
        % integrals anyways...
        function[val] = delta_function_x(obj)
            val = obj.mag_grad_phi_cd .* obj.delta_function_phi;
        end
        
        % return a bump function on the object's grid
        function[val] = bump_function(obj, x0, y0, rad)
            val = LevelSet.bump_function_2d(obj.x, obj.y, x0, y0, rad);
        end
        
        % "cd" = central differencing
        function[dphi_dx, dphi_dy] = grad_phi_cd(obj)
            dphi_dx = 0.5 * (obj.phi_pad(4:end-3,5:end-2) - obj.phi_pad(4:end-3,3:end-4));
            dphi_dy = 0.5 * (obj.phi_pad(5:end-2,4:end-3) - obj.phi_pad(3:end-4,4:end-3));
        end
        
        function[val] = mag_grad_phi_cd(obj)
            [dphi_dx, dphi_dy] = grad_phi_cd(obj);
            val = sqrt(dphi_dx.*dphi_dx + dphi_dy.*dphi_dy);
        end
        
        % for normals, compute one-sided differences at the
        % boundaries, since the normals are often not continuous (the
        % level-set function is not smooth) at a periodic boundary
        function[nx, ny] = normals(obj, x, y)
            [dphi_dx, dphi_dy] = obj.grad_phi_cd();
            dphi_dx(:,1)   = obj.phi_pad(4:end-3,5)     - obj.phi_pad(4:end-3,4);
            dphi_dx(:,end) = obj.phi_pad(4:end-3,end-3) - obj.phi_pad(4:end-3,end-4);
            dphi_dy(1,:)   = obj.phi_pad(5,4:end-3)     - obj.phi_pad(4,4:end-3);
            dphi_dy(end,:) = obj.phi_pad(end-3,4:end-3) - obj.phi_pad(end-4,4:end-3);
            nx = dphi_dx ./ (sqrt(dphi_dx.*dphi_dx + dphi_dy.*dphi_dy) + 1e-10);
            ny = dphi_dy ./ (sqrt(dphi_dx.*dphi_dx + dphi_dy.*dphi_dy) + 1e-10);
            if (nargin>1)
                Fx = griddedInterpolant(obj.x.', obj.y.', nx.');
                nx = Fx(x(:), y(:));
                Fy = griddedInterpolant(obj.x.', obj.y.', ny.');
                ny = Fy(x(:), y(:));
            end
        end
        
        function[val] = delta_function_phi(obj)
            val = 0.5 / obj.smooth_eps * (1 + cos(pi * obj.phi / obj.smooth_eps)) .* (abs(obj.phi) <= obj.smooth_eps);
        end
        
        function[H] = step_function_phi(obj)
            H = 0.5 * (1 + obj.phi / obj.smooth_eps + sin(pi * obj.phi / obj.smooth_eps) / pi) .* (abs(obj.phi) <= obj.smooth_eps) ...
                + 1 * (obj.phi > obj.smooth_eps);
        end
        
        % This derivative of this smoothed Heaviside function gives the
        % delta function above
        function[val] = int_vol_integral(obj, fn_vol)
            val = sum(sum(fn_vol .* (1 - obj.step_function_phi)));
        end
        
        function[val] = ext_vol_integral(obj, fn_vol)
            H = 0.5 * (1 + obj.phi / obj.smooth_eps + sin(pi * obj.phi / obj.smooth_eps) / pi) .* (abs(obj.phi) <= obj.smooth_eps) ...
                + 1 * (obj.phi > obj.smooth_eps);
            val = sum(sum(fn_vol .* H));
        end
        
        function[val] = vol(obj)
            val = obj.int_vol_integral(ones(size(obj.phi)));
        end
        
        function[val] = SA(obj)
            val = obj.surface_integral(ones(size(obj.phi)));
        end
        
        % phi_pad holds phi in a larger matrix with periodic padding,
        %   to enable faster indexing when taking differences (without
        %   "for" loops)
        function[obj] = update_phi_pad(obj)
            Nx = obj.g.Nx;
            Ny = obj.g.Ny;
            pad_sz = 3; % for 5th-order ENO/WENO
            obj.phi_pad = zeros(Ny+2*pad_sz, Nx+2*pad_sz);
            obj.phi_pad(1+pad_sz:Ny+pad_sz, 1+pad_sz:Nx+pad_sz) = obj.phi;
            obj.phi_pad(1+pad_sz:Ny+pad_sz, 1:pad_sz) = obj.phi(:,4:-1:2);
            obj.phi_pad(1+pad_sz:Ny+pad_sz, Nx+1+pad_sz:Nx+2*pad_sz) = obj.phi(:,end-1:-1:end-3);
            obj.phi_pad(1:pad_sz,:) = obj.phi_pad(7:-1:5,:);
            obj.phi_pad(Ny+1+pad_sz:Ny+2*pad_sz,:) = obj.phi_pad(Ny+pad_sz-1:-1:Ny+pad_sz-3,:); 
        end
        
    end
    
    methods (Static)
        % Sec. 3.4 in Fedkiw & Osher
        function[val] = bump_function_2d(x, y, x0, y0, rad)
            val = LevelSet.bump_function_1d(x,x0,rad) .* LevelSet.bump_function_1d(y,y0,rad);
        end
        
        function[val] = bump_function_1d(x, x0, rad)
            val = zeros(size(x));
            cond = (abs((x-x0)/rad)<1);
            val(cond) = exp(-1 ./ (1 - abs((x(cond) - x0)/rad).^2));
%             val = val;
        end
        
        function dphi = hj_weno_dphi(v1, v2, v3, v4, v5)
            S1 = 13 / 12 * (v1 - 2*v2 + v3).^2 + 0.25 * (v1 - 4*v2 + 3*v3).^2;
            S2 = 13 / 12 * (v2 - 2*v3 + v4).^2 + 0.25 * (v2 - v4).^2;
            S3 = 13 / 12 * (v3 - 2*v4 + v5).^2 + 0.25 * (3*v3 - 4*v4 + v5).^2;
            
            epsilon = 1e-6 * max([v1(:)'.^2; v2(:)'.^2; v3(:)'.^2; v4(:)'.^2; v5(:)'.^2])' + 1e-99;
            epsilon = reshape(epsilon', size(v1));
            alpha1 = 0.1 ./ (S1 + epsilon).^2;
            alpha2 = 0.6 ./ (S2 + epsilon).^2;
            alpha3 = 0.3 ./ (S3 + epsilon).^2;
            w1 = alpha1 ./ (alpha1 + alpha2 + alpha3);
            w2 = alpha2 ./ (alpha1 + alpha2 + alpha3);
            w3 = alpha3 ./ (alpha1 + alpha2 + alpha3);
            
            dphi = w1 .* (v1/3 - 7*v2/6 + 11*v3/6) ...
                  + w2 .* (-v2/6 + 5*v3/6 + v4/3) ...
                  + w3 .* (v3/3 + 5*v4/6 - v5/6);            
        end
        
        function[sd] = circle(x, y, x0, y0, r)
            sd = sqrt((x-x0).^2 + (y-y0).^2) - r;
        end
        
        function[sd] = rectangle(x, y, x0, y0, x1, y1)
            v0 = [x1 - x0; 0];
            v1 = [0; y1 - y0];
            xy_vec = [x(:) - 0.5*(x0+x1), y(:) - 0.5*(y0+y1)];
            a = abs(xy_vec * v0) / dot(v0, v0);
            b = abs(xy_vec * v1) / dot(v1, v1);
            a0 = min(a,0.5);
            b0 = min(b,0.5);
            sd = -1 * min((0.5 - a) * norm(v0), (0.5 - b) * norm(v1)) .* logical((a<0.5).*(b<0.5)) ...
                + sqrt((a-a0).^2.*dot(v0,v0) + (b-b0).^2.*dot(v1,v1));
        end
        
%         % ax, ay are lengths of unit cell in x, y directions
%         % by symmetry no need for periodic circle (other shapes would need
%         %   periodic repetition)
%         function[sd] = periodic_circle(x, y, x0, y0, r, ax, ay)
%             sd0 = sqrt((x-x0).^2 + (y-y0).^2) - r;
%             sd1 = min(sd0, sqrt((x-x0-ax).^2 + (y-y0).^2) - r);
%             sd2 = min(sd1, sqrt((x-x0+ax).^2 + (y-y0).^2) - r);
%             sd3 = min(sd2, sqrt((x-x0).^2 + (y-y0-ay).^2) - r);
%             sd = min(sd3, sqrt((x-x0).^2 + (y-y0+ay).^2) - r);
%         end
    end
end

