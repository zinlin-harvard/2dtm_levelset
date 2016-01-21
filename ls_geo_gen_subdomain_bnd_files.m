
addpath('~/Level-Set-Opt/mesh');

alpha = 0.8; % courant number
phi = @(x,y) LevelSet.circle(x, y, 0.5, 0.5, 0.4);

Nx = 200;
Ny = Nx;
lx = 1; % sizes of unit cell
ly = 1; % note that xN - x1 = lx - dx, similarly for y

dx = lx / Nx;
dy = ly / Ny;
grid = grid2d(Nx, Ny, 0, 0, dx, dy);
grido = offset_grid(grid);

ls = LevelSet(grido, alpha, phi);

bc_xy = {'pec', 'per'};
pml_xy = [1 1];
lpml_xy = [0.3 0.2];
lsg = LSGeo(ls, bc_xy, pml_xy, lpml_xy);

lsg.gen_mesh('mesh.msh', 'constraints.txt');
m = Mesh('mesh.msh');
m = m.reorder_triangles();
domains = lsg.mark_domains_xml(m, 'subdomains.xml');
lsg.mark_bnds_xml(m, 'bnds.xml', domains);

GRID=[ls.g.x,ls.g.y];
save('grid.txt','GRID','-ascii');
[nx,ny] = lsg.ls.normals();
nx=transpose(nx);
ny=transpose(ny);
nx=nx(:);
ny=ny(:);
phidata=transpose(ls.phi);
phidata=phidata(:);
dlmwrite('grid_ls_normals.txt',[ls.g.x ls.g.y nx ny], 'delimiter', ' ', 'precision', 16);
save('grid_ls_phi.txt','phidata','-ascii');