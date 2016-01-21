function Vn=compute_velocity(lambda,eps1,eps2,Nx,Ny)

Eraw=load('EzfieldGrid.txt');
phi=load('grid_ls_phi.txt');
omega=2*pi/lambda;

phishift=(1-sign(phi))/2;
epsilon=(eps2-eps1)*phishift+eps1;

Ez=Eraw(:,3)+ 1i*Eraw(:,4);
Ezadj=Ez/(1i*omega);

Vn = real( (eps2-eps1) * Ez.*Ezadj );

end