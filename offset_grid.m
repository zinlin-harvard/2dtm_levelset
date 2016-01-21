
% dx and dx are optional (defaults are 0.5*dx, 0.5*dy, as for FDQS offset)
function[go] = offset_grid(g, deltax, deltay)
    if(nargin<3)
        deltax = 0.5 * g.dx;
        deltay = 0.5 * g.dy;
    end
    go = grid2d(g.Nx, g.Ny, g.x0 + deltax, g.y0 + deltay, g.dx, g.dy);
end