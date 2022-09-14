%% Data
posmay=5;
EL = zz50.mix_energy_list;
DF = zz50.dens_fermi;
DB = zz50.dens_bec;


%% plotting the deltaE in Fermi
deltae = zeros(2553,1);
cc=0;
for j=1:3000
    cc = cc+1;
    piv = EL(:,j);
    deltae(cc) = piv(Ntg+1)-piv(Ntg);    
end

figure
plot(deltae,'.')
xlabel('g_m'); ylabel('\Delta E');

title("Gap Energy Evolution in maximum position of y ="+posmay)
saveas(gcf,"Horizontal_phase_"+Ntg+"_particles_y_"+posmay+".png")
saveas(gcf,"Horizontal_phase_"+Ntg+"_particles_y_"+posmay+".fig")
%% plotting the gm evolution
cc = 0;
figure
for ii=1:2500
    m = mod(ii,500);
    if ii==1 || m==0 %|| ii==2253 %|| ii==2254
        cc=cc+1;
        subplot(3,2,cc)
        plot(EL(1:Ntg+5,ii),'o')
        xlabel('Number of particles'); ylabel('Energy');
        title("g_m = " + round(ii/1000,2))
    end
end

saveas(gcf,"Horizontal_gp_"+Ntg+"_particles_y_"+posmay+".png")
%% Density Fermi Particles Plot
cc = 0;

x1 = linspace(-posmax,posmax,Nx);
y1 = linspace(-posmay,posmay,Ny);
[x2d, y2d] = meshgrid(x1,y1);

figure
for ii=1:2500
    m = mod(ii,500);
    if ii==1 || m==0 %|| ii==2250 %|| ii==2254
        cc=cc+1;
        subplot(3,2,cc)
        pcolor(x2d, y2d, DF(:,:,ii));
        xlabel('x'); ylabel('y ');
        subtitle("Fermi Gas in g_m = " + round(ii/1000,2))
        shading flat; axis square; colorbar;
    end
end
saveas(gcf,"Horizontal_fp_"+Ntg+"_particles_y_"+posmay+".png")

%% Density Fermi Video
cc = 0;

x1 = linspace(-posmax,posmax,Nx);
y1 = linspace(-posmay,posmay,Ny);
[x2d, y2d] = meshgrid(x1,y1);

pcolor(x2d,y2d,DF(:,:,1));
xlabel('x'); ylabel('y');
subtitle("g_m = " + j/1000)
shading flat; axis square; colorbar;
fr=getframe;


myVideo = VideoWriter("evolution_fdens_"+Ntg+"_particles_y_"+posmay); %open video file
myVideo.FrameRate = 5;  %can adjust this, 5 - 10 works well for me

open(myVideo)
for k=1:3000
    if k==1 || mod(k,100)==0 
        var = DF(:,:,k);
        pcolor(x2d,y2d,var);
        xlabel('x'); ylabel('y');
        subtitle("Fermi Gas in g_m = " + k/1000)
        shading flat; axis square; colorbar;
        %fr(k) = getframe;
        drawnow
        pause(0.4)
        frame = getframe(gcf); %get 
        writeVideo(myVideo, frame);
    end
end
close(myVideo)

%% Density BEC Particles Plot
cc = 0;

x1 = linspace(-posmax,posmax,Nx);
y1 = linspace(-posmay,posmay,Ny);
[x2d, y2d] = meshgrid(x1,y1);

figure
for ii=1:2500
    m = mod(ii,500);
    if ii==1 || m==0 %|| ii==2250 %|| ii==2254
        cc=cc+1;
        subplot(3,2,cc)
        pcolor(x2d, y2d, DB(:,:,ii));
        xlabel('x'); ylabel('y ');
        subtitle("BEC in g_m = " + round(ii/1000,2))
        shading flat; axis square; colorbar;
    end
end
saveas(gcf,"Horizontal_bp_"+Ntg+"_particles_y_"+posmay+".png")

%% Density BEC Video
cc = 0;

x1 = linspace(-posmax,posmax,Nx);
y1 = linspace(-posmay,posmay,Ny);
[x2d, y2d] = meshgrid(x1,y1);


myVideo = VideoWriter("evolution_bdens_"+Ntg+"_particles_y_"+posmay); %open video file
myVideo.FrameRate = 5;  %can adjust this, 5 - 10 works well for me

open(myVideo)
for k=1:3000
    if k==1 || mod(k,100)==0 
        var = DB(:,:,k);
        pcolor(x2d,y2d,var);
        xlabel('x'); ylabel('y');
        subtitle("BEC in g_m = " + k/1000)
        shading flat; axis square; colorbar;
        %fr(k) = getframe;
        drawnow
        pause(0.4)
        frame = getframe(gcf); %get 
        writeVideo(myVideo, frame);
    end
end
close(myVideo)
