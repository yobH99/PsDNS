clc; clear all; close all;
caso = "Re_1000_16_32"; 
folder_data_spectra = '/home/yobh/Desktop/LANL/PsDNS/caso_buoyant/'+caso+'/spectra_2D';

addpath(folder_data_spectra)


filenames = import_names(folder_data_spectra) ;
nfiles = length(filenames);
[nk,nz,nt,zvec,tvec] = extract_parameter(filenames);

Eu  = zeros(nk,nz,nt);
Euk = zeros(nk,nz,nt);
Ev  = zeros(nk,nz,nt);
Evk = zeros(nk,nz,nt);
Ew  = zeros(nk,nz,nt);
Ewk = zeros(nk,nz,nt);
Ec  = zeros(nk,nz,nt);
Eck = zeros(nk,nz,nt);


for n = 2:nfiles
    tmp = importdata(filenames{n});
    kvec = tmp(:,1);

    Eu  (:,zvec(n)+1,tvec(n)+1) = tmp(:,2);
    Euk (:,zvec(n)+1,tvec(n)+1) = tmp(:,3);
    Ev  (:,zvec(n)+1,tvec(n)+1) = tmp(:,4);
    Evk (:,zvec(n)+1,tvec(n)+1) = tmp(:,5);
    Ew  (:,zvec(n)+1,tvec(n)+1) = tmp(:,6);
    Ewk (:,zvec(n)+1,tvec(n)+1) = tmp(:,7);
    Ec  (:,zvec(n)+1,tvec(n)+1) = tmp(:,8);
    Eck (:,zvec(n)+1,tvec(n)+1) = tmp(:,9);
    
end
%%
fig_folder = "/home/yobh/Desktop/UofM/AEROSP_525/project/figures"
close all
%time instants
nt = length(Ec(1,1,1:30));

gray_levels = linspace(0.9, 0.01, nt);
colors = repmat(gray_levels', 1, 3);
%colors(:,3) = 0.2
%colors(:,1) = 0.5

sizes = [500, 500, 800, 600]; 
ftsize = 24;


fig = figure('Position', sizes, 'Visible', 'on');
ax = axes(fig);

for t = 5:nt
    loglog(ax, kvec, abs(Ec(:,44,t)) .* kvec.^(5/3), ...
        'LineWidth', 1.5,'color',colors(t,:))
    hold on 
end
xlabel(ax, '$k$', 'Interpreter', 'latex', 'FontSize', ftsize);
ylabel(ax, '$k^{5/3} E_c(k,0)$', 'Interpreter', 'latex', 'FontSize', ftsize);
set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', ftsize);
grid(ax, 'on')
box(ax, 'on')

print(fig, fig_folder+'/spectrum_plot_'+caso+'.eps', '-depsc')
%close(fig)

%%
clc
clear all
close all

caso = "Re_1000_32_64";
folder_data = '/home/yobh/Desktop/LANL/PsDNS/caso_buoyant/' + caso + '/data';
addpath(folder_data)

tic
data_dissipation = import_profiles('dissipation.dat');
data_profiles    = import_profiles('profiles.dat');
nz = 355;

[nz,z, epsxx, epsyy, epszz, epsxy, epsxz, S, G, ...
 u, v, w, c, ...
 Rxx, Ryy, Rzz, Rxy, Rxz, Ryz, ...
 cc, ax, ay, az, ...
 Rxxx, Ryyy, Rzzz, Rxxy, Rxxz, Ryyx, Ryyz, Rzzx, Rzzy, Rxyz, ...
 Rcxx, Rcyy, Rczz, Rcxy, Rcxz, Rcyz, Rccx, Rccy, Rccz] = load_all_data(data_dissipation, data_profiles);

Lz = 8*pi;
zvec = linspace(-Lz, Lz, nz);


%%
close all 
clc
z0 = 20;
dt = 0.5;
[nz,nt] = size(c);
nt = floor(nt/2);
tvec = (1:nt)*dt;
ctmp = (c(z0:end-z0,:)+1)*0.5;
ztmp = zvec(z0:end-z0);
eps = 0.01;

hs = zeros(nt,1);
hb = zeros(nt,1);

for t = 1:nt
    %bubble height 

    idxs = find(ctmp(:,t) > eps, 1, 'first')
    idxb = find(ctmp(:,t) < 1-eps,1,'last');
    ctmp(:,t)
    hs(t) = zvec(idxs);
    hb(t) = zvec(idxb);
 
end
%hb = hb - hb(1);
%hs = hs - hs(1);
plot(tvec,hs,tvec,hb)
vis = 0.003;
Reb = fd4(hb,dt).*hb/vis;
Res = fd4(hs,dt).*hs/vis;


h = hb - hs% 


figure(2)
loglog(tvec,h,'*-k')
hold on 
loglog(tvec(floor(nt/2)+1:end), tvec(floor(nt/2)+1:end).^2/10)
hold on 
loglog(tvec(5:floor(nt/3)+1), 10*tvec(5:floor(nt/3)+1).^0.5/10)





function [nz, z, epsxx, epsyy, epszz, epsxy, epsxz, S, G, ...
          u, v, w, c, ...
          Rxx, Ryy, Rzz, Rxy, Rxz, Ryz, ...
          cc, ax, ay, az, ...
          Rxxx, Ryyy, Rzzz, Rxxy, Rxxz, Ryyx, Ryyz, Rzzx, Rzzy, Rxyz, ...
          Rcxx, Rcyy, Rczz, Rcxy, Rcxz, Rcyz, Rccx, Rccy, Rccz] = ...
          load_all_data(data_dissipation, data_profiles)
    
    nz = 355;  % If you know it, otherwise extract automatically
    nblocks = size(data_dissipation, 1) / nz;
    
    if mod(nblocks,1) ~= 0
        error('Data size mismatch: total rows not divisible by nz.');
    end

    nblocks = round(nblocks);

    z     = zeros(nz, nblocks);
    epsxx = zeros(nz, nblocks);
    epsyy = zeros(nz, nblocks);
    epszz = zeros(nz, nblocks);
    epsxy = zeros(nz, nblocks);
    epsxz = zeros(nz, nblocks);
    S     = zeros(nz, nblocks);
    G     = zeros(nz, nblocks);
    
    u     = zeros(nz, nblocks);
    v     = zeros(nz, nblocks);
    w     = zeros(nz, nblocks);
    c     = zeros(nz, nblocks);
    
    Rxx   = zeros(nz, nblocks);
    Ryy   = zeros(nz, nblocks);
    Rzz   = zeros(nz, nblocks);
    Rxy   = zeros(nz, nblocks);
    Rxz   = zeros(nz, nblocks);
    Ryz   = zeros(nz, nblocks);
    
    cc    = zeros(nz, nblocks);
    ax    = zeros(nz, nblocks);
    ay    = zeros(nz, nblocks);
    az    = zeros(nz, nblocks);
    
    Rxxx  = zeros(nz, nblocks);
    Ryyy  = zeros(nz, nblocks);
    Rzzz  = zeros(nz, nblocks);
    Rxxy  = zeros(nz, nblocks);
    Rxxz  = zeros(nz, nblocks);
    Ryyx  = zeros(nz, nblocks);
    Ryyz  = zeros(nz, nblocks);
    Rzzx  = zeros(nz, nblocks);
    Rzzy  = zeros(nz, nblocks);
    Rxyz  = zeros(nz, nblocks);
    
    Rcxx  = zeros(nz, nblocks);
    Rcyy  = zeros(nz, nblocks);
    Rczz  = zeros(nz, nblocks);
    Rcxy  = zeros(nz, nblocks);
    Rcxz  = zeros(nz, nblocks);
    Rcyz  = zeros(nz, nblocks);
    Rccx  = zeros(nz, nblocks);
    Rccy  = zeros(nz, nblocks);
    Rccz  = zeros(nz, nblocks);
    
    startidx = 1;

    for t = 1:nblocks
        endidx = startidx + nz - 1;
    
        % If needed you can uncomment this
        % z(:,t) = data_dissipation(startidx:endidx,1);

        epsxx(:,t) = data_dissipation(startidx:endidx,2);
        epsyy(:,t) = data_dissipation(startidx:endidx,3);
        epszz(:,t) = data_dissipation(startidx:endidx,4);
        epsxy(:,t) = data_dissipation(startidx:endidx,5);
        epsxz(:,t) = data_dissipation(startidx:endidx,6);
        S(:,t)     = data_dissipation(startidx:endidx,7);
        G(:,t)     = data_dissipation(startidx:endidx,8);
    
        u(:,t)     = data_profiles(startidx:endidx,2);
        v(:,t)     = data_profiles(startidx:endidx,3);
        w(:,t)     = data_profiles(startidx:endidx,4);
        c(:,t)     = data_profiles(startidx:endidx,5);
    
        Rxx(:,t)   = data_profiles(startidx:endidx,6);
        Ryy(:,t)   = data_profiles(startidx:endidx,7);
        Rzz(:,t)   = data_profiles(startidx:endidx,8);
        Rxy(:,t)   = data_profiles(startidx:endidx,9);
        Rxz(:,t)   = data_profiles(startidx:endidx,10);
        Ryz(:,t)   = data_profiles(startidx:endidx,11);
    
        cc(:,t)    = data_profiles(startidx:endidx,12);
        ax(:,t)    = data_profiles(startidx:endidx,13);
        ay(:,t)    = data_profiles(startidx:endidx,14);
        az(:,t)    = data_profiles(startidx:endidx,15);
    
        Rxxx(:,t)  = data_profiles(startidx:endidx,16);
        Ryyy(:,t)  = data_profiles(startidx:endidx,17);
        Rzzz(:,t)  = data_profiles(startidx:endidx,18);
        Rxxy(:,t)  = data_profiles(startidx:endidx,19);
        Rxxz(:,t)  = data_profiles(startidx:endidx,20);
        Ryyx(:,t)  = data_profiles(startidx:endidx,21);
        Ryyz(:,t)  = data_profiles(startidx:endidx,22);
        Rzzx(:,t)  = data_profiles(startidx:endidx,23);
        Rzzy(:,t)  = data_profiles(startidx:endidx,24);
        Rxyz(:,t)  = data_profiles(startidx:endidx,25);
    
        Rcxx(:,t)  = data_profiles(startidx:endidx,26);
        Rcyy(:,t)  = data_profiles(startidx:endidx,27);
        Rczz(:,t)  = data_profiles(startidx:endidx,28);
        Rcxy(:,t)  = data_profiles(startidx:endidx,29);
        Rcxz(:,t)  = data_profiles(startidx:endidx,30);
        Rcyz(:,t)  = data_profiles(startidx:endidx,31);
        Rccx(:,t)  = data_profiles(startidx:endidx,32);
        Rccy(:,t)  = data_profiles(startidx:endidx,33);
        Rccz(:,t)  = data_profiles(startidx:endidx,34);
    
        startidx = endidx + 1;
    end
end

function filenames = import_names(folder_data)
      if exist(folder_data, 'dir')
          files = dir(fullfile(folder_data, '*'));
          files = files(~ismember({files.name}, {'.', '..'}));
          filenames = {files.name};
          disp('Files found in the folder:');
      else
          disp('Folder does not exist.');
      end
end

function [zNum,tNum] = extract_zt(filename)
      zPattern = 'z(\d+)';
      tPattern = 't(\d+)';
      zTokens = regexp(filename, zPattern, 'tokens');
      tTokens = regexp(filename, tPattern, 'tokens');

      if ~isempty(zTokens)
          zNum = str2double(zTokens{1}{1});
      else
          zNum = NaN;
      end

      if ~isempty(tTokens)
          tToken = tTokens{1}{1};
          extractedT = '';
          encounteredFirstZero = false;
          for i = 1:length(tToken)
              if ~encounteredFirstZero && tToken(i)=='0'
                  encounteredFirstZero = true;
                  extractedT = [extractedT tToken(i)];
              elseif encounteredFirstZero && tToken(i)=='0'
                  break;
              else
                  extractedT = [extractedT tToken(i)];
              end
          end
          tNum = str2double(extractedT);
      else
          tNum = NaN;
      end
end

function [nk,nz,nt,zvec,tvec] = extract_parameter(filenames)
      nfiles = length(filenames);
      nk = max(size(importdata(filenames{2})));
      tvec = zeros(1,nfiles);
      zvec = zeros(1,nfiles);
      for n = 2:nfiles
          filename = filenames{n};
          [zNum, tNum] = extract_zt(filename);
          zvec(n) = int32(zNum);
          tvec(n) = int32(tNum);
      end
      nz = max(zvec);
      nt = max(tvec);
end

function dudx = fd4(u, dx)

n = length(u);
dudx = zeros(size(u));

dudx(1) = (-25*u(1) + 48*u(2) - 36*u(3) + 16*u(4) - 3*u(5)) / (12*dx);
dudx(2) = (-3*u(1) - 10*u(2) + 18*u(3) - 6*u(4) + u(5)) / (12*dx);

for i = 3:n-2
    dudx(i) = (u(i-2) - 8*u(i-1) + 8*u(i+1) - u(i+2)) / (12*dx);
end

dudx(n-1) = (3*u(n) + 10*u(n-1) - 18*u(n-2) + 6*u(n-3) - u(n-4)) / (12*dx);
dudx(n) = (25*u(n) - 48*u(n-1) + 36*u(n-2) - 16*u(n-3) + 3*u(n-4)) / (12*dx);

end