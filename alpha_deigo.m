%% Mix of a Fermionic Gas in a BEC in adiabatic evolution
%%Generic code for an arbitrary potential


%%--------BEC and Fermionic gas parameters

nx = 5;                                 %number of computational basis in x axis
ny = 5;                                 %number of computational basis in y axis
Nx = 2^nx;                              % x-grid points number
Ny = 2^ny;                              % y-grid points number
posmax =5;                              %maximum number in the grid space in x --> in 15 is the position where all particles are counted
posmaylist = linspace(0.1,1,10);%1;     %maximum number in the grid space in y--> in 15 is the position where all particles are counted
posmay = posmaylist(zz);
dt = 1e-5;                              %step size 
NBEC=1e4;                               %Number of atoms
Ntg =3;                                 %Number of fermions
steps = 50;                             %time step for the evolution
omegax = 1;                             %harmonic trap in x --> if working with box potential turn it zero
omegay=1;                               %harmonic trap in y --> if working with box potential turn it zero

[x,dx,px,dpx]=fftdef(posmax,Nx);        % Creating the grid points in x position
[y,dy,py,dpy]=fftdef(posmay,Ny);        % Creating the grid points in y position
 
[ym,xm]=meshgrid(y,x);                  %meshgrid for x and y     
r=sqrt(xm.^2 + ym.^2);                  %position vector in cylyndrical coordinates
[pym,pxm]=meshgrid(py,px);              %meshgrid for px and py 
p=sqrt(pxm.^2 + pym.^2);                %momentum vector in cylyndrical coordinates

gc=1;                                   % Intraspecies interaction term
maxgm = 2;                              %maximum number of points
pgm = 1000;                             %points per unit
gm=linspace(0,maxgm,maxgm*pgm);         % Interspecies interaction term


%%--------------BEC Hamiltonian operators 

%Vbec=0.5*(bomegax^2.*xm.^2 + bomegay^2.*ym.^2);        %Harmonic Trap
Vbec = 0;                                               %Thomas-Fermi approximation
Kbec = exp(-0.5*dt*p.^2);                               %Kinetic Term

%%------------- Fermionic Gas Operators
Ekinx = -1/(2*dx^2)*(diag(ones(Nx-1,1),-1) - 2*diag(ones(Nx,1)) + diag(ones(Nx-1,1),1));        % approximating the kinetic energy term with a finite difference scheme in x 
Ekiny = -1/(2*dy^2)*(diag(ones(Ny-1,1),-1) - 2*diag(ones(Ny,1)) + diag(ones(Ny-1,1),1));        % approximating the kinetic energy term with a finite difference scheme in y
%Vtrapx = diag(0.5*(fomegax^2)*x.^2);                                                           % harmonic trapping potential
%Vtrapy = diag(0.5*(fomegay^2)*y.^2);

idx = diag(ones(Nx,1));                                                                         %identity matrix                
idy = diag(ones(Ny,1));                                                                         %identity matrix
kinxy = kron(Ekinx,idy)+kron(idx,Ekiny);                                                        %Kinetic Potential in function of x and y --> H_x(x)I_x+I_y(x)H_y 

%Vtrapxy= 0.5*((fomegax^2).*(xm.^2) +(fomegay^2).*(ym.^2));                                     %Potential Part in function of x and y  
%Vpotxy = diag(reshape(Vtrapxy,[],1));                                                          %reshaping in n^2 times n^2 dimmension
Vpotxy=0;
fham=kinxy+Vpotxy;                                                                              %fermionic hamiltonian

%%----------------------- Calculating the Initial Wave function-------------------------

[ES,EV] = eig(fham);                                                                            %diagonalizing the Hamiltonian 
energies = diag(EV);                                                                            %extracting the eigenenergies
states = (1/sqrt(dx*dy)).*ES;                                                                   %normalizing the eigenstates

initial_fermi_energy = sum(energies(1:Ntg));                                                    %sum of the lowest Ntg fermi particles energies

init_ss = sum(abs(states(:,1:Ntg)).^2,2);                                                       %summing the lowest Ntg wave  functions
denswft=reshape(init_ss,Ny,Nx);                                                                 %density of the wavepacket of the lowest wavefunctions of the gas

%a0x=6;                                                                                         %calibration of the harmonic trap in x direction
%a0y=6;                                                                                         %calibration of the harmonic trap in y direction
%wfc = real(((mu - Vbec)/(gc)).^(1/2));                                                         %Thomas-Fermi wave function
%wfc=exp(-((xm./(2*a0x)).^2 + (ym./(2*a0y)).^2));                                               %Gaussian wavefunction guess.

wfc=5*ones(Ny,Nx);                                                                              %Flat wavefunction
wfc=sqrt(NBEC).*wfc./(sum(sum(abs(wfc).^2))*dx*dy)^(0.5);                                       %Renormalise for groundstate

mix_energy = 0;                                                                                 %Initializing our space where energy will be saved    
count = 0;

%%-----------------------------Imaginary time-------------------------------
plotcount=0;
%figure
gm_list=zeros(maxgm*pgm,1);
mix_energy_list=zeros(Nx*Ny,maxgm*pgm);
dens_fermi=zeros(Ny,Nx,maxgm*pgm);
%bec_energy_list=zeros(maxgm*pgm,1);
%fermi_energy_list=zeros(maxgm*pgm,1);


%figure
for ii=gm
    for aa=1:1000  
        for becloop=1:steps
        
            V = Vbec + gc*abs(wfc).^2 + ii*denswft;                             %Creating the Potential term, sum of the trapping and the intraspecies interaction                                           
            
            wfc = exp(-0.5*dt*V).*wfc;                                          %Combining the wave function with a half of the potential part
            
            wfc = ifft2(Kbec.'.*fft2(wfc));                                     %Fast Fourier Transform
            
            wfc = exp(-0.5*dt*V).*wfc;                                          %Combining the wave function with a half of the potential part
            
            wfc=(sqrt(NBEC)).*wfc./(sum(sum(abs(wfc).^2))*dx*dy).^(0.5);        %Renormalise for groundstate
            
            %wfb = wfc;                                                         %Keeping the value in another variable to further plotting                               
        end
        %%---- Adding the fermionic Gas
    
        dens_bec = abs(wfc).^2;
        inter=diag(reshape(ii*dens_bec,[],1));
    
    
        [ES,EV] = eig(fham+inter);                                              %diagonalizing the Hamiltonian 
        energies = diag(EV);                                                    %extracting the eigenenergies
        states = (1/sqrt(dx*dy)).*ES;                                           %normalizing the eigenstates
    
    
        %%--------Taking the fermionic density of the lowest Ntg waves
                                                                                %Number of Ntg particles
        ss = sum(abs(states(:,1:Ntg)).^2,2);                                    %summing the lowest Ntg wave  functions
        denswft=reshape(ss,Ny,Nx);                                              %density of the wavepacket of the lowest wavefunctions of the gas
    
        % Energy of the Ntg low particles                                    
        fermi_energy = sum(energies(1:Ntg));                                    %sum of the lowest Ntg fermi particles energies
    
        % Energies for BEC
        Ken = 0.5*conj(wfc).*ifft2(((p.').^2).*fft2(wfc));
        Ven = Vbec.*abs(wfc).^2 + 0.5*gc.*abs(wfc).^4;
        bec_energy = sum(sum(abs(Ken + Ven))*dx)*dy;
    
        %Total Energy
        total_energy = bec_energy+fermi_energy;
    
        count=count+1;
    
        if abs(total_energy - mix_energy) < 1e-3                                %Finding a truncted ground state
            total_energy;
            break
        end
        mix_energy=total_energy;
                  
    end
    plotcount = plotcount +1;
    
    mix_energy_list(:,plotcount) = energies;
    dens_fermi(:,:,plotcount)=denswft;
    %plotcount
end

filepath = 'Temporary/';
filename=strcat('zz=',string(zz),'.mat');

save(strcat(filepath,filename))