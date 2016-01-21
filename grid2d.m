classdef grid2d
    %grid2d contains mapping from 1:N to points on grid, neighbors
    %   Detailed explanation goes here
    
    properties
        % values
        x0 
        y0
        dx
        dy
        Nx
        Ny
        N
        
        % arrays
        x 
        y
        
        % indices
        rb % right boundary
        lb % left
        bb % top
        tb % bottom
    end
    
    methods
        function obj = grid2d(Nx, Ny, x0, y0, dx, dy)       
            obj.Nx = Nx;
            obj.Ny = Ny;
            obj.N = Nx * Ny;
            obj.x0 = x0;
            obj.y0 = y0;
            obj.dx = dx;
            obj.dy = dy;

            obj.x = repmat((1:Nx).', [1 Ny]);
            obj.y = repmat((1:Ny), [Nx 1]);
            
            obj.x = obj.dx * (obj.x(:) - 1) + x0;
            obj.y = obj.dy * (obj.y(:) - 1) + y0;
                        
            j = Nx;
            i = (1:Ny).';
            obj.rb = j + (i-1)*Nx;
            
            i = Ny;
            j = (1:Nx).';
            obj.tb = j + (i-1)*Nx;
            
            j = 1;
            i = (1:Ny).';
            obj.lb = j + (i-1)*Nx;
            
            i = 1;
            j = (1:Nx).';
            obj.bb = j + (i-1)*Nx;
        end 
        
        % this function assumes within one period of unit cell
        function[ind] = per_ind(obj,j,i)
            i = i + obj.Ny.*(i<1) - obj.Ny.*(i>obj.Ny);
            j = j + obj.Nx.*(j<1) - obj.Nx.*(j>obj.Nx);
            
            ind = j + (i-1)*obj.Nx;
        end
    end
end

% original code for per_ind
%             if (i<1)
%                 i = i + obj.Ny;
%             elseif (i>obj.Ny)
%                 i = i - obj.Ny;
%             end
%             if (j<1)
%                 j = j + obj.Nx;
%             elseif (j>obj.Nx)
%                 j = j - obj.Nx;
%             end

