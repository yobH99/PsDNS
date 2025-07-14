clc 
clear all 
close all 
Lz = 4*pi;
nz = 355;
cmap = [ ...
    0.20, 0.00, 0.30;   % deep purple
    1.00, 0.20, 0.60;   % hot pink
    1.00, 0.60, 0.10;   % neon orange
    0.20, 0.40, 1.00;   % electric blue
    0.50, 0.75, 1.00    % sky cyan
];

set(groot,'DefaultAxesColorOrder',    cmap)


%basic settings
sizes = [500, 500, 800, 600];
ftsize = 32;
vis    = 0.0035;

% list of cases
%cases     = { "Re_1000_0_8", "Re_1000_16_32", "Re_1000_32_64" };
cases = {"Re_1000_40_60_pi_pi_4pi"};
num_cases = numel(cases);

% storage for mixing lengths
h_all   = cell(num_cases,1);
t_all   = cell(num_cases,1);
hb_all  = cell(num_cases,1);
hs_all  = cell(num_cases,1);
Re_all  = cell(num_cases,1);

% storage for correlation lengths
lambdaU_all = cell(num_cases,1);
lambdaV_all = cell(num_cases,1);
lambdaW_all = cell(num_cases,1);
lambdaC_all = cell(num_cases,1);

% paths (adjust if needed)
base_spectra = "/home/yobh/Desktop/LANL/PsDNS/caso_buoyant/";
base_data    = "/home/yobh/Desktop/LANL/PsDNS/caso_buoyant/";
fig_folder   = "/home/yobh/Desktop/UofM/AEROSP_525/project/figures/";

%loop over cases
for i = 1:num_cases
    caso = cases{i};
    
    % 1) load and assemble spectra
    folder_spec = fullfile(base_spectra, caso, 'spectra_2D');
    addpath(folder_spec)
    filenames   = import_names(folder_spec);
    [~,~,nt,zvec,tvec_spec] = extract_parameter(filenames);
    nt = floor(nt/3);
    zphys=linspace(-Lz/2,Lz/2,floor(nz/4)+1);

    % clear any old arrays
    clear Eu Euk Ev Evk Ew Ewk Ec Eck
    
    for n = 2:numel(filenames)
        tmp = importdata(fullfile(folder_spec, filenames{n}));
        kvec = tmp(:,1);
        Eu  (:, zvec(n)+1, tvec_spec(n)+1) = abs(tmp(:,2));
        Euk (:, zvec(n)+1, tvec_spec(n)+1) = abs(tmp(:,3));
        Ev  (:, zvec(n)+1, tvec_spec(n)+1) = abs(tmp(:,4));
        Evk (:, zvec(n)+1, tvec_spec(n)+1) = abs(tmp(:,5));
        Ew  (:, zvec(n)+1, tvec_spec(n)+1) = abs(tmp(:,6));
        Ewk (:, zvec(n)+1, tvec_spec(n)+1) = abs(tmp(:,7));
        Ec  (:, zvec(n)+1, tvec_spec(n)+1) = abs(tmp(:,8));
        Eck(:, zvec(n)+1, tvec_spec(n)+1) = abs(tmp(:,9));
    end

    % pick midpoint slice
    midz = round(mean(zvec)) + 1;
    
    % compute correlation lengths
    lambdaU = lambda(kvec,zphys,Eck,Ec);
    lambdaV = lambda(kvec,zphys, Evk, Ev);
    lambdaW = lambda(kvec,zphys, Ewk, Ew);
    lambdaC = lambda(kvec,zphys, Eck, Ec);
    
    lambdaU_all{i} = lambdaU;
    lambdaV_all{i} = lambdaV;
    lambdaW_all{i} = lambdaW;
    lambdaC_all{i} = lambdaC;
    
    % 1b) plot compensated concentration spectra
    figSpec = figure('Position', sizes, 'Visible', 'on');
    axSpec  = axes(figSpec);
    colors  = inferno(nt);
    for t = 40:nt
        loglog(axSpec, kvec, Ec(:,midz,t).*kvec.^(7/4), ...
               'LineWidth',1.5, 'Color', colors(t,:),'LineStyle','-');
        hold(axSpec,'on');
    end
    hold(axSpec,'off');
    xlabel(axSpec,'$k$','Interpreter','latex','FontSize',ftsize);
    ylabel(axSpec,'$k^{7/4}E_c(k,0)$','Interpreter','latex','FontSize',ftsize);
    set(axSpec,'TickLabelInterpreter','latex','FontSize',ftsize);
    ylim([1e-11,inf])
    grid(axSpec,'on'); box(axSpec,'on');
    print(figSpec, fullfile(fig_folder, "spectrum_" + caso + ".eps"), '-depsc');
    
    % 2) mixing‐layer thickness from profiles
    folder_dat = fullfile(base_data, caso, 'data');
    addpath(folder_dat);

    data_d = import_profiles('dissipation.dat');
    data_p = import_profiles('profiles.dat');
    [nz, z, epsxx, epsyy, epszz, epsxy, epsxz, S, G, ...
     u, v, w, c, Rxx, Ryy, Rzz, Rxy, Rxz, Ryz, ...
     cc, ax, ay, az, Rxxx, Ryyy, Rzzz, Rxxy, Rxxz, ...
     Ryyx, Ryyz, Rzzx, Rzzy, Rxyz, Rcxx, Rcyy, Rczz, ...
     Rcxy, Rcxz, Rcyz, Rccx, Rccy, Rccz] = load_all_data(data_d, data_p);

    nz   = size(c,1);
    dt   = 0.5;
    z0   = 20;
    nt_m = floor(size(c,2)/3);

    zvec_full = linspace(-8*pi, 8*pi, nz);
    ctmp = (c(z0:end-z0, :) + 1) * 0.5;
    ztmp = zvec_full(z0:end-z0);

    eps = 0.01;
    hs  = zeros(nt_m,1);
    hb  = zeros(nt_m,1);

    for t = 1:nt_m
        idx_s = find(ctmp(:,t) >  eps,    1, 'first');
        idx_b = find(ctmp(:,t) < 1 - eps, 1, 'last');
        hs(t) = ztmp(idx_s);
        hb(t) = ztmp(idx_b);
    end

    tvec_h    = (1:nt_m)' * dt;
    hb_all{i} = hb;
    hs_all{i} = hs;
    h_all{i}  = hb - hs;
    t_all{i}  = tvec_h;
end

%% overlay mixing-layer thickness from all cases
fig1 = figure('Position', sizes, 'Visible', 'on');
ax1  = axes(fig1);
hold(ax1,'on')
legend_names = {'$k_m \in [0,8]$', '$k_m \in [16,32]$', '$k_m \in [32,64]$'};
for i = 1:num_cases
    plot(ax1, t_all{i}, h_all{i}, 'LineWidth', 2);
end
hold(ax1,'off')
legend(ax1, legend_names, 'Interpreter','latex','Location','best')
xlabel(ax1, '$t$', 'Interpreter','latex','FontSize',ftsize)
ylabel(ax1, 'mixing length $h(t)$', 'Interpreter','latex','FontSize',24)
set(ax1, 'TickLabelInterpreter','latex','FontSize',ftsize)
grid(ax1,'on'); box(ax1,'on')
print(fig1, fullfile(fig_folder, "mixing_length.eps"), '-depsc')

%% Re_h over time
fig2 = figure('Position', sizes, 'Visible', 'on');
ax2  = axes(fig2);
hold(ax2,'on')
for i = 1:num_cases
    Re_all{i} = abs(h_all{i} .* fd4(h_all{i}, dt) / vis);
    plot(ax2, t_all{i}, Re_all{i}, 'LineWidth', 2);
end
hold(ax2,'off')
legend(ax2, legend_names, 'Interpreter','latex','Location','best')
xlabel(ax2, '$t$', 'Interpreter','latex','FontSize',ftsize)
ylabel(ax2, '$\mathrm{Re}_h(t)$', 'Interpreter','latex','FontSize',ftsize)
set(ax2, 'TickLabelInterpreter','latex','FontSize',ftsize)
grid(ax2,'on'); box(ax2,'on')
xlim(ax2, [0,8])
print(fig2, fullfile(fig_folder, 'mixing_length_Reh.eps'), '-depsc')

%% 4) one figure per λ–type across all cases (use S1, S2, S3)
lambda_cell   = { lambdaU_all,    lambdaV_all,    lambdaW_all,    lambdaC_all };
lambda_labels = { 'U',            'V',            'W',            'C' };
case_labels   = { 'S1',           'S2',           'S3' };
dt  = 0.5;
for j = 1:numel(lambda_cell)
    fig = figure('Position', sizes, 'Visible', 'on');
    ax  = axes(fig);
    %hold(ax,'on');
    
    for i = 1:num_cases
        t_idx = 1:numel(lambda_cell{j}{i});
        loglog(ax, t_idx*dt, lambda_cell{j}{i}, 'LineWidth',1.5, ...
               'DisplayName', sprintf('$\\lambda_{%s}$, %s', ...
                          lambda_labels{j}, case_labels{i}));
        hold on 
    end
    %ttmp = t_idx(15:30)*dt/2;
    %loglog(ax, ttmp,ttmp.^0.5)
    
    hold(ax,'off');
    legend(ax, 'Interpreter','latex', 'Location','northwest');
    xlabel(ax,'$t$', 'Interpreter','latex','FontSize',ftsize);
    ylabel(ax, sprintf('$\\lambda_{%s}(t)$', lambda_labels{j}), ...
           'Interpreter','latex','FontSize',ftsize);
    set(ax,'TickLabelInterpreter','latex','FontSize',ftsize);
    grid(ax,'on'); box(ax,'on');
    
    print(fig, fullfile(fig_folder, ...
          sprintf('lambda_%s_allcases.eps', lower(lambda_labels{j}))), '-depsc');
end



function corr_length = lambda(kvec,zvec,Eck,Ec)   
    corr_length_z =  squeeze(trapz(kvec,Eck, 1)./trapz(kvec,Ec,1));
    size(corr_length_z)
    corr_length   =  trapz(zvec,corr_length_z,1);
end


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


function [cm_data]=inferno(m)
cm = [[  1.46159096e-03,   4.66127766e-04,   1.38655200e-02],
       [  2.26726368e-03,   1.26992553e-03,   1.85703520e-02],
       [  3.29899092e-03,   2.24934863e-03,   2.42390508e-02],
       [  4.54690615e-03,   3.39180156e-03,   3.09092475e-02],
       [  6.00552565e-03,   4.69194561e-03,   3.85578980e-02],
       [  7.67578856e-03,   6.13611626e-03,   4.68360336e-02],
       [  9.56051094e-03,   7.71344131e-03,   5.51430756e-02],
       [  1.16634769e-02,   9.41675403e-03,   6.34598080e-02],
       [  1.39950388e-02,   1.12247138e-02,   7.18616890e-02],
       [  1.65605595e-02,   1.31362262e-02,   8.02817951e-02],
       [  1.93732295e-02,   1.51325789e-02,   8.87668094e-02],
       [  2.24468865e-02,   1.71991484e-02,   9.73274383e-02],
       [  2.57927373e-02,   1.93306298e-02,   1.05929835e-01],
       [  2.94324251e-02,   2.15030771e-02,   1.14621328e-01],
       [  3.33852235e-02,   2.37024271e-02,   1.23397286e-01],
       [  3.76684211e-02,   2.59207864e-02,   1.32232108e-01],
       [  4.22525554e-02,   2.81385015e-02,   1.41140519e-01],
       [  4.69146287e-02,   3.03236129e-02,   1.50163867e-01],
       [  5.16437624e-02,   3.24736172e-02,   1.59254277e-01],
       [  5.64491009e-02,   3.45691867e-02,   1.68413539e-01],
       [  6.13397200e-02,   3.65900213e-02,   1.77642172e-01],
       [  6.63312620e-02,   3.85036268e-02,   1.86961588e-01],
       [  7.14289181e-02,   4.02939095e-02,   1.96353558e-01],
       [  7.66367560e-02,   4.19053329e-02,   2.05798788e-01],
       [  8.19620773e-02,   4.33278666e-02,   2.15289113e-01],
       [  8.74113897e-02,   4.45561662e-02,   2.24813479e-01],
       [  9.29901526e-02,   4.55829503e-02,   2.34357604e-01],
       [  9.87024972e-02,   4.64018731e-02,   2.43903700e-01],
       [  1.04550936e-01,   4.70080541e-02,   2.53430300e-01],
       [  1.10536084e-01,   4.73986708e-02,   2.62912235e-01],
       [  1.16656423e-01,   4.75735920e-02,   2.72320803e-01],
       [  1.22908126e-01,   4.75360183e-02,   2.81624170e-01],
       [  1.29284984e-01,   4.72930838e-02,   2.90788012e-01],
       [  1.35778450e-01,   4.68563678e-02,   2.99776404e-01],
       [  1.42377819e-01,   4.62422566e-02,   3.08552910e-01],
       [  1.49072957e-01,   4.54676444e-02,   3.17085139e-01],
       [  1.55849711e-01,   4.45588056e-02,   3.25338414e-01],
       [  1.62688939e-01,   4.35542881e-02,   3.33276678e-01],
       [  1.69575148e-01,   4.24893149e-02,   3.40874188e-01],
       [  1.76493202e-01,   4.14017089e-02,   3.48110606e-01],
       [  1.83428775e-01,   4.03288858e-02,   3.54971391e-01],
       [  1.90367453e-01,   3.93088888e-02,   3.61446945e-01],
       [  1.97297425e-01,   3.84001825e-02,   3.67534629e-01],
       [  2.04209298e-01,   3.76322609e-02,   3.73237557e-01],
       [  2.11095463e-01,   3.70296488e-02,   3.78563264e-01],
       [  2.17948648e-01,   3.66146049e-02,   3.83522415e-01],
       [  2.24762908e-01,   3.64049901e-02,   3.88128944e-01],
       [  2.31538148e-01,   3.64052511e-02,   3.92400150e-01],
       [  2.38272961e-01,   3.66209949e-02,   3.96353388e-01],
       [  2.44966911e-01,   3.70545017e-02,   4.00006615e-01],
       [  2.51620354e-01,   3.77052832e-02,   4.03377897e-01],
       [  2.58234265e-01,   3.85706153e-02,   4.06485031e-01],
       [  2.64809649e-01,   3.96468666e-02,   4.09345373e-01],
       [  2.71346664e-01,   4.09215821e-02,   4.11976086e-01],
       [  2.77849829e-01,   4.23528741e-02,   4.14392106e-01],
       [  2.84321318e-01,   4.39325787e-02,   4.16607861e-01],
       [  2.90763373e-01,   4.56437598e-02,   4.18636756e-01],
       [  2.97178251e-01,   4.74700293e-02,   4.20491164e-01],
       [  3.03568182e-01,   4.93958927e-02,   4.22182449e-01],
       [  3.09935342e-01,   5.14069729e-02,   4.23720999e-01],
       [  3.16281835e-01,   5.34901321e-02,   4.25116277e-01],
       [  3.22609671e-01,   5.56335178e-02,   4.26376869e-01],
       [  3.28920763e-01,   5.78265505e-02,   4.27510546e-01],
       [  3.35216916e-01,   6.00598734e-02,   4.28524320e-01],
       [  3.41499828e-01,   6.23252772e-02,   4.29424503e-01],
       [  3.47771086e-01,   6.46156100e-02,   4.30216765e-01],
       [  3.54032169e-01,   6.69246832e-02,   4.30906186e-01],
       [  3.60284449e-01,   6.92471753e-02,   4.31497309e-01],
       [  3.66529195e-01,   7.15785403e-02,   4.31994185e-01],
       [  3.72767575e-01,   7.39149211e-02,   4.32400419e-01],
       [  3.79000659e-01,   7.62530701e-02,   4.32719214e-01],
       [  3.85228383e-01,   7.85914864e-02,   4.32954973e-01],
       [  3.91452659e-01,   8.09267058e-02,   4.33108763e-01],
       [  3.97674379e-01,   8.32568129e-02,   4.33182647e-01],
       [  4.03894278e-01,   8.55803445e-02,   4.33178526e-01],
       [  4.10113015e-01,   8.78961593e-02,   4.33098056e-01],
       [  4.16331169e-01,   9.02033992e-02,   4.32942678e-01],
       [  4.22549249e-01,   9.25014543e-02,   4.32713635e-01],
       [  4.28767696e-01,   9.47899342e-02,   4.32411996e-01],
       [  4.34986885e-01,   9.70686417e-02,   4.32038673e-01],
       [  4.41207124e-01,   9.93375510e-02,   4.31594438e-01],
       [  4.47428382e-01,   1.01597079e-01,   4.31080497e-01],
       [  4.53650614e-01,   1.03847716e-01,   4.30497898e-01],
       [  4.59874623e-01,   1.06089165e-01,   4.29845789e-01],
       [  4.66100494e-01,   1.08321923e-01,   4.29124507e-01],
       [  4.72328255e-01,   1.10546584e-01,   4.28334320e-01],
       [  4.78557889e-01,   1.12763831e-01,   4.27475431e-01],
       [  4.84789325e-01,   1.14974430e-01,   4.26547991e-01],
       [  4.91022448e-01,   1.17179219e-01,   4.25552106e-01],
       [  4.97257069e-01,   1.19379132e-01,   4.24487908e-01],
       [  5.03492698e-01,   1.21575414e-01,   4.23356110e-01],
       [  5.09729541e-01,   1.23768654e-01,   4.22155676e-01],
       [  5.15967304e-01,   1.25959947e-01,   4.20886594e-01],
       [  5.22205646e-01,   1.28150439e-01,   4.19548848e-01],
       [  5.28444192e-01,   1.30341324e-01,   4.18142411e-01],
       [  5.34682523e-01,   1.32533845e-01,   4.16667258e-01],
       [  5.40920186e-01,   1.34729286e-01,   4.15123366e-01],
       [  5.47156706e-01,   1.36928959e-01,   4.13510662e-01],
       [  5.53391649e-01,   1.39134147e-01,   4.11828882e-01],
       [  5.59624442e-01,   1.41346265e-01,   4.10078028e-01],
       [  5.65854477e-01,   1.43566769e-01,   4.08258132e-01],
       [  5.72081108e-01,   1.45797150e-01,   4.06369246e-01],
       [  5.78303656e-01,   1.48038934e-01,   4.04411444e-01],
       [  5.84521407e-01,   1.50293679e-01,   4.02384829e-01],
       [  5.90733615e-01,   1.52562977e-01,   4.00289528e-01],
       [  5.96939751e-01,   1.54848232e-01,   3.98124897e-01],
       [  6.03138930e-01,   1.57151161e-01,   3.95891308e-01],
       [  6.09330184e-01,   1.59473549e-01,   3.93589349e-01],
       [  6.15512627e-01,   1.61817111e-01,   3.91219295e-01],
       [  6.21685340e-01,   1.64183582e-01,   3.88781456e-01],
       [  6.27847374e-01,   1.66574724e-01,   3.86276180e-01],
       [  6.33997746e-01,   1.68992314e-01,   3.83703854e-01],
       [  6.40135447e-01,   1.71438150e-01,   3.81064906e-01],
       [  6.46259648e-01,   1.73913876e-01,   3.78358969e-01],
       [  6.52369348e-01,   1.76421271e-01,   3.75586209e-01],
       [  6.58463166e-01,   1.78962399e-01,   3.72748214e-01],
       [  6.64539964e-01,   1.81539111e-01,   3.69845599e-01],
       [  6.70598572e-01,   1.84153268e-01,   3.66879025e-01],
       [  6.76637795e-01,   1.86806728e-01,   3.63849195e-01],
       [  6.82656407e-01,   1.89501352e-01,   3.60756856e-01],
       [  6.88653158e-01,   1.92238994e-01,   3.57602797e-01],
       [  6.94626769e-01,   1.95021500e-01,   3.54387853e-01],
       [  7.00575937e-01,   1.97850703e-01,   3.51112900e-01],
       [  7.06499709e-01,   2.00728196e-01,   3.47776863e-01],
       [  7.12396345e-01,   2.03656029e-01,   3.44382594e-01],
       [  7.18264447e-01,   2.06635993e-01,   3.40931208e-01],
       [  7.24102613e-01,   2.09669834e-01,   3.37423766e-01],
       [  7.29909422e-01,   2.12759270e-01,   3.33861367e-01],
       [  7.35683432e-01,   2.15905976e-01,   3.30245147e-01],
       [  7.41423185e-01,   2.19111589e-01,   3.26576275e-01],
       [  7.47127207e-01,   2.22377697e-01,   3.22855952e-01],
       [  7.52794009e-01,   2.25705837e-01,   3.19085410e-01],
       [  7.58422090e-01,   2.29097492e-01,   3.15265910e-01],
       [  7.64009940e-01,   2.32554083e-01,   3.11398734e-01],
       [  7.69556038e-01,   2.36076967e-01,   3.07485188e-01],
       [  7.75058888e-01,   2.39667435e-01,   3.03526312e-01],
       [  7.80517023e-01,   2.43326720e-01,   2.99522665e-01],
       [  7.85928794e-01,   2.47055968e-01,   2.95476756e-01],
       [  7.91292674e-01,   2.50856232e-01,   2.91389943e-01],
       [  7.96607144e-01,   2.54728485e-01,   2.87263585e-01],
       [  8.01870689e-01,   2.58673610e-01,   2.83099033e-01],
       [  8.07081807e-01,   2.62692401e-01,   2.78897629e-01],
       [  8.12239008e-01,   2.66785558e-01,   2.74660698e-01],
       [  8.17340818e-01,   2.70953688e-01,   2.70389545e-01],
       [  8.22385784e-01,   2.75197300e-01,   2.66085445e-01],
       [  8.27372474e-01,   2.79516805e-01,   2.61749643e-01],
       [  8.32299481e-01,   2.83912516e-01,   2.57383341e-01],
       [  8.37165425e-01,   2.88384647e-01,   2.52987700e-01],
       [  8.41968959e-01,   2.92933312e-01,   2.48563825e-01],
       [  8.46708768e-01,   2.97558528e-01,   2.44112767e-01],
       [  8.51383572e-01,   3.02260213e-01,   2.39635512e-01],
       [  8.55992130e-01,   3.07038188e-01,   2.35132978e-01],
       [  8.60533241e-01,   3.11892183e-01,   2.30606009e-01],
       [  8.65005747e-01,   3.16821833e-01,   2.26055368e-01],
       [  8.69408534e-01,   3.21826685e-01,   2.21481734e-01],
       [  8.73740530e-01,   3.26906201e-01,   2.16885699e-01],
       [  8.78000715e-01,   3.32059760e-01,   2.12267762e-01],
       [  8.82188112e-01,   3.37286663e-01,   2.07628326e-01],
       [  8.86301795e-01,   3.42586137e-01,   2.02967696e-01],
       [  8.90340885e-01,   3.47957340e-01,   1.98286080e-01],
       [  8.94304553e-01,   3.53399363e-01,   1.93583583e-01],
       [  8.98192017e-01,   3.58911240e-01,   1.88860212e-01],
       [  9.02002544e-01,   3.64491949e-01,   1.84115876e-01],
       [  9.05735448e-01,   3.70140419e-01,   1.79350388e-01],
       [  9.09390090e-01,   3.75855533e-01,   1.74563472e-01],
       [  9.12965874e-01,   3.81636138e-01,   1.69754764e-01],
       [  9.16462251e-01,   3.87481044e-01,   1.64923826e-01],
       [  9.19878710e-01,   3.93389034e-01,   1.60070152e-01],
       [  9.23214783e-01,   3.99358867e-01,   1.55193185e-01],
       [  9.26470039e-01,   4.05389282e-01,   1.50292329e-01],
       [  9.29644083e-01,   4.11479007e-01,   1.45366973e-01],
       [  9.32736555e-01,   4.17626756e-01,   1.40416519e-01],
       [  9.35747126e-01,   4.23831237e-01,   1.35440416e-01],
       [  9.38675494e-01,   4.30091162e-01,   1.30438175e-01],
       [  9.41521384e-01,   4.36405243e-01,   1.25409440e-01],
       [  9.44284543e-01,   4.42772199e-01,   1.20354038e-01],
       [  9.46964741e-01,   4.49190757e-01,   1.15272059e-01],
       [  9.49561766e-01,   4.55659658e-01,   1.10163947e-01],
       [  9.52075421e-01,   4.62177656e-01,   1.05030614e-01],
       [  9.54505523e-01,   4.68743522e-01,   9.98735931e-02],
       [  9.56851903e-01,   4.75356048e-01,   9.46952268e-02],
       [  9.59114397e-01,   4.82014044e-01,   8.94989073e-02],
       [  9.61292850e-01,   4.88716345e-01,   8.42893891e-02],
       [  9.63387110e-01,   4.95461806e-01,   7.90731907e-02],
       [  9.65397031e-01,   5.02249309e-01,   7.38591143e-02],
       [  9.67322465e-01,   5.09077761e-01,   6.86589199e-02],
       [  9.69163264e-01,   5.15946092e-01,   6.34881971e-02],
       [  9.70919277e-01,   5.22853259e-01,   5.83674890e-02],
       [  9.72590351e-01,   5.29798246e-01,   5.33237243e-02],
       [  9.74176327e-01,   5.36780059e-01,   4.83920090e-02],
       [  9.75677038e-01,   5.43797733e-01,   4.36177922e-02],
       [  9.77092313e-01,   5.50850323e-01,   3.90500131e-02],
       [  9.78421971e-01,   5.57936911e-01,   3.49306227e-02],
       [  9.79665824e-01,   5.65056600e-01,   3.14091591e-02],
       [  9.80823673e-01,   5.72208516e-01,   2.85075931e-02],
       [  9.81895311e-01,   5.79391803e-01,   2.62497353e-02],
       [  9.82880522e-01,   5.86605627e-01,   2.46613416e-02],
       [  9.83779081e-01,   5.93849168e-01,   2.37702263e-02],
       [  9.84590755e-01,   6.01121626e-01,   2.36063833e-02],
       [  9.85315301e-01,   6.08422211e-01,   2.42021174e-02],
       [  9.85952471e-01,   6.15750147e-01,   2.55921853e-02],
       [  9.86502013e-01,   6.23104667e-01,   2.78139496e-02],
       [  9.86963670e-01,   6.30485011e-01,   3.09075459e-02],
       [  9.87337182e-01,   6.37890424e-01,   3.49160639e-02],
       [  9.87622296e-01,   6.45320152e-01,   3.98857472e-02],
       [  9.87818759e-01,   6.52773439e-01,   4.55808037e-02],
       [  9.87926330e-01,   6.60249526e-01,   5.17503867e-02],
       [  9.87944783e-01,   6.67747641e-01,   5.83286889e-02],
       [  9.87873910e-01,   6.75267000e-01,   6.52570167e-02],
       [  9.87713535e-01,   6.82806802e-01,   7.24892330e-02],
       [  9.87463516e-01,   6.90366218e-01,   7.99897176e-02],
       [  9.87123759e-01,   6.97944391e-01,   8.77314215e-02],
       [  9.86694229e-01,   7.05540424e-01,   9.56941797e-02],
       [  9.86174970e-01,   7.13153375e-01,   1.03863324e-01],
       [  9.85565739e-01,   7.20782460e-01,   1.12228756e-01],
       [  9.84865203e-01,   7.28427497e-01,   1.20784651e-01],
       [  9.84075129e-01,   7.36086521e-01,   1.29526579e-01],
       [  9.83195992e-01,   7.43758326e-01,   1.38453063e-01],
       [  9.82228463e-01,   7.51441596e-01,   1.47564573e-01],
       [  9.81173457e-01,   7.59134892e-01,   1.56863224e-01],
       [  9.80032178e-01,   7.66836624e-01,   1.66352544e-01],
       [  9.78806183e-01,   7.74545028e-01,   1.76037298e-01],
       [  9.77497453e-01,   7.82258138e-01,   1.85923357e-01],
       [  9.76108474e-01,   7.89973753e-01,   1.96017589e-01],
       [  9.74637842e-01,   7.97691563e-01,   2.06331925e-01],
       [  9.73087939e-01,   8.05409333e-01,   2.16876839e-01],
       [  9.71467822e-01,   8.13121725e-01,   2.27658046e-01],
       [  9.69783146e-01,   8.20825143e-01,   2.38685942e-01],
       [  9.68040817e-01,   8.28515491e-01,   2.49971582e-01],
       [  9.66242589e-01,   8.36190976e-01,   2.61533898e-01],
       [  9.64393924e-01,   8.43848069e-01,   2.73391112e-01],
       [  9.62516656e-01,   8.51476340e-01,   2.85545675e-01],
       [  9.60625545e-01,   8.59068716e-01,   2.98010219e-01],
       [  9.58720088e-01,   8.66624355e-01,   3.10820466e-01],
       [  9.56834075e-01,   8.74128569e-01,   3.23973947e-01],
       [  9.54997177e-01,   8.81568926e-01,   3.37475479e-01],
       [  9.53215092e-01,   8.88942277e-01,   3.51368713e-01],
       [  9.51546225e-01,   8.96225909e-01,   3.65627005e-01],
       [  9.50018481e-01,   9.03409063e-01,   3.80271225e-01],
       [  9.48683391e-01,   9.10472964e-01,   3.95289169e-01],
       [  9.47594362e-01,   9.17399053e-01,   4.10665194e-01],
       [  9.46809163e-01,   9.24168246e-01,   4.26373236e-01],
       [  9.46391536e-01,   9.30760752e-01,   4.42367495e-01],
       [  9.46402951e-01,   9.37158971e-01,   4.58591507e-01],
       [  9.46902568e-01,   9.43347775e-01,   4.74969778e-01],
       [  9.47936825e-01,   9.49317522e-01,   4.91426053e-01],
       [  9.49544830e-01,   9.55062900e-01,   5.07859649e-01],
       [  9.51740304e-01,   9.60586693e-01,   5.24203026e-01],
       [  9.54529281e-01,   9.65895868e-01,   5.40360752e-01],
       [  9.57896053e-01,   9.71003330e-01,   5.56275090e-01],
       [  9.61812020e-01,   9.75924241e-01,   5.71925382e-01],
       [  9.66248822e-01,   9.80678193e-01,   5.87205773e-01],
       [  9.71161622e-01,   9.85282161e-01,   6.02154330e-01],
       [  9.76510983e-01,   9.89753437e-01,   6.16760413e-01],
       [  9.82257307e-01,   9.94108844e-01,   6.31017009e-01],
       [  9.88362068e-01,   9.98364143e-01,   6.44924005e-01]];
if nargin < 1
    cm_data = cm;
else
    hsv=rgb2hsv(cm);
    hsv(144:end,1)=hsv(144:end,1)+1; % hardcoded
    cm_data=interp1(linspace(0,1,size(cm,1)),hsv,linspace(0,1,m));
    cm_data(cm_data(:,1)>1,1)=cm_data(cm_data(:,1)>1,1)-1;
    cm_data=hsv2rgb(cm_data);
  
end
end