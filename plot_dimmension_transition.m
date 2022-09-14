%% Plotting the results
zz1 = load("zz=1.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz2 = load("zz=2.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz3 = load("zz=3.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz4 = load("zz=4.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz5 = load("zz=5.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz6 = load("zz=6.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz7 = load("zz=7.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz8 = load("zz=8.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz9 = load("zz=9.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz10 = load("zz=10.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz11 = load("zz=11.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz12 = load("zz=12.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz13 = load("zz=13.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz14 = load("zz=14.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz15 = load("zz=15.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz16 = load("zz=16.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz17 = load("zz=17.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz18 = load("zz=18.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz19 = load("zz=19.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz20 = load("zz=20.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz21 = load("zz=21.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz22 = load("zz=22.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz23 = load("zz=23.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz24 = load("zz=24.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz25 = load("zz=25.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz26 = load("zz=26.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz27 = load("zz=27.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz28 = load("zz=28.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz29 = load("zz=29.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz30 = load("zz=30.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz31 = load("zz=31.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz32 = load("zz=32.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz33 = load("zz=33.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz34 = load("zz=34.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz35 = load("zz=35.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz36 = load("zz=36.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz37 = load("zz=37.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz38 = load("zz=38.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz39 = load("zz=39.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz40 = load("zz=40.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz41 = load("zz=41.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz42 = load("zz=42.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz43 = load("zz=43.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz44 = load("zz=44.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz45 = load("zz=45.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz46 = load("zz=46.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz47 = load("zz=47.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz48 = load("zz=48.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz49 = load("zz=49.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');
zz50 = load("zz=50.mat",'mix_energy_list', 'dens_fermi', 'dens_bec');

z =[zz1,zz2,zz3,zz4,zz5,zz6,zz7,zz8,zz9,zz10...
    ,zz11,zz12,zz13,zz14,zz15,zz16,zz17,zz18,zz19,zz20...
    ,zz21,zz22,zz23,zz24,zz25,zz26,zz27,zz28,zz29,zz30...
    ,zz31,zz32,zz33,zz34,zz35,zz36,zz37,zz38,zz39,zz40...
    ,zz41,zz42,zz43,zz44,zz45,zz46,zz47,zz48,zz49,zz50];
%%
list_energy = zeros(1024,3000,50);
list_densfermi = zeros(32,32,3000,50);
list_densbec = zeros(32,32,3000,50);

for i=1:50
    list_energy(:,:,i)=z(i).mix_energy_list;
    list_densfermi(:,:,:,i)=z(i).dens_fermi;
    list_densbec(:,:,:,i)=z(i).dens_bec;
end

%% Variables

Ntg = 7;
posmax=5;
            %Variables from 1 to 3000
%% plotting the deltaE
gm=1;          
deltae = zeros(50,1);
xxx=linspace(1,50,50);
cc=0;
for j=10:50
    cc = cc+1;
    piv = list_energy(:,gm,j); %(1:Ntg+10,gm,ii)
    deltae(cc) = piv(Ntg+1)-piv(Ntg);   
end

figure
plot(xxx./10,deltae,'o')
xlabel('posmay'); ylabel('\Delta E');
title("Gap Energy Evolution in variable y-maximum position at g_m = "+ gm/1000)
saveas(gcf,"Gap_energy_gm_"+gm/1000+"_particles_"+Ntg+".png")
saveas(gcf,"Gap_energy_gm_"+gm/1000+"_particles_"+Ntg+".fig")



%% Density of the BEC
gm=2300;
cc = 0;
Nx=32;
Ny=32;

XB = list_densbec;
myVideo = VideoWriter("bec_density at g_m_"+ gm/1000+"_particles_"+Ntg); %open video file
myVideo.FrameRate = 10;  %can adjust this, 5 - 10 works well for me

open(myVideo)
for k=22:50
    if mod(k,1)==0 
        posmay = k/10;
        x1 = linspace(-posmax,posmax,Nx);
        y1 = linspace(-posmay,posmay,Ny);
        [x2d, y2d] = meshgrid(x1,y1);
        var = XB(:,:,gm,k);
        pcolor(x2d,y2d,var);
        xlabel('x'); ylabel('y');
        subtitle("posmay = " + k/10);
        shading flat; axis square; colorbar;
        %fr(k) = getframe;
        drawnow
        pause(0.4)
        frame = getframe(gcf); %get 
        writeVideo(myVideo, frame);
    end
end
close(myVideo)


%% Density of the Fermi
gm=2300;
cc = 0;
Nx=32;
Ny=32;

XX = list_densfermi;
myVideo = VideoWriter("fermi_density at g_m_"+ gm/1000+"_particles_"+Ntg); %open video file
myVideo.FrameRate = 10;  %can adjust this, 5 - 10 works well for me

open(myVideo)
for k=22:50
    if mod(k,1)==0 
        posmay = k/10;
        x1 = linspace(-posmax,posmax,Nx);
        y1 = linspace(-posmay,posmay,Ny);
        [x2d, y2d] = meshgrid(x1,y1);
        var = XX(:,:,gm,k);
        pcolor(x2d,y2d,var);
        xlabel('x'); ylabel('y');
        subtitle("posmay = " + k/10);
        shading flat; axis square; colorbar;
        %fr(k) = getframe;
        drawnow
        pause(0.4)
        frame = getframe(gcf); %get 
        writeVideo(myVideo, frame);
    end
end
close(myVideo)

%% 

%list_energy(:,:,i)
X = list_energy;
Ntg=3;
h=plot(X(1:Ntg+10,gm,1),'o');

xlabel('Number of particles'); ylabel('Energy');
%title("posmay = " );
fr=getframe; %I don't think there will be much benefit to a more careful initialization here...

%%Initialize video
myVideo = VideoWriter("energies at gm_"+ (gm/1000)); %open video file
myVideo.FrameRate = 5;  %can adjust this, 5 - 10 works well for me

open(myVideo)
for ii = 1:50
    if mod(ii,1)==0 
        h.YData = X(1:Ntg+10,gm,ii);
        %fr(i) = getframe; 
        % or if you just want to watch it and you don't care about the playback
        % speed, consider 'drawnow limitrate'
        xlabel('Number of particles'); ylabel('Energy');
        title("posmay = " + ii/10);
        pause(0.25)
        frame = getframe(gcf); %get 
        writeVideo(myVideo, frame);
    end   
end
close(myVideo)
% write a video file, or call movie





