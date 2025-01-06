
% *********************************************************************
% MATLAB Script for Calculating Bootstrap Current Coefficients
% *********************************************************************
%
% This MATLAB script computes the bootstrap current coefficients based 
% on input data from the 'neo_input.nc' file. 
%
% Output Files:
% --------------------
% The script generates one of the following output file names:
% 
% 1. ProfileJBSCoeff_Te_L31_32_34_alpha_B2_dtedpsit_G
% 2. ProfileJBSCoeff_Psi_L31_32_34_alpha_B2_dtedpsit_G
% 
% The output files contain the following 8 columns:
% 
% Column 1: Te or Psi         
% Column 2: L31               - Bootstrap current coefficient L31.
% Column 3: L32               - Bootstrap current coefficient L32.
% Column 4: L34               - Bootstrap current coefficient L34.
% Column 5: alpha             - Amplification factor for the bootstrap current.
% Column 6: 1/<B²>            - Inverse of Flux surface average magnetic field squared.
% Column 7: dTe/dPsi_t * 2π   - 
% Column 8: Gbar/(i-N)        - Landreman et al (2022)'s isomorphism parameter
%
% These coefficients are calculated offline using the 'neo_input.nc' file as input.
% The script also plots these coefficients, uncomment this part below to
% generate the plots
% *********************************************************************
%
% File Placement:
% -----------------
% The output file must be placed in the M3D-C1 working case directory 
% for further processing in M3D-C1 simulations.
%
% *********************************************************************
%
% M3D-C1 Input File Settings:
% ------------------------------
% In order to use the bootstrap current coefficients in M3D-C1, the 
% following parameters must be set in the M3D-C1 input file:
%
% - ibootstrap = 2            ! 1 when using profiles as a function of psi
%                             ! 2 when using dTe/dPsi_t for bootstrap 
% - ibootstrap_model = 4      ! Choose the bootstrap model:
%                             !  1 = Sauter (with eqsubtract = 0/1)
%                             !  2 = Redl   (with eqsubtract = 0/1)
%                             !  3 = Sauter (with eqsubtract = 0)
%                             !  4 = Redl   (with eqsubtract = 0)
% - bootstrap_alpha = 1       ! Amplification factor for the bootstrap current.
%                             !  Set to 1 for no amplification.
% - ibootstrap_map_te = 1     ! 1 to map bootstrap coefficients with 
%                             !  respect to Te (used for stellarator cases).
%                             !  This setting is to be used with 
%                             !  ibootstrap = 2, ibootstrap_model = 3/4.
% - ibootstrap_regular = 1e-8 ! default 1e-8
%                             ! regularization term for the expression = ((del_p, del_Te)) / (grad_Te_magnitude**2 + adaptive_regularization)

% *********************************************************************

% Main script
clc;
clear;
close all;



%%
%========================================================
% Set initial flags and parameters
%========================================================
%
% - UseSauter0Redl1 = 1       ! Choose bootstrap current model:
%                             !  0 = Sauter Calculations
%                             !  1 = Redl Calculations
% - imap_te = 1               !  1: mapping of the bootstrap coefficients w.r.t Te,
% - N = 0                     !  N from Landreman's formula: 
%                             !   Gbar/(i-N) = (G + NI) / (i - N), where N is the number of grid points.
% - UseBfore = 1              % 1: epsilon = (Bmax - Bmin) / (Bmax + Bmin);
%                                 0: epsilon = Rmin/Rc %Rc =(R_+ + R_-)/2 see IIIA of Hager & Chang (2016)
%                                 where R_+ and R_- are major radii at the outer and inner-most points of a flux surface         
% - psi_start, psi_end    
% - filename = neo_input.nc filename
% *********************************************************************

UseSauter0Redl1 = 1; %0: Sauter Calculations, 2: Redl Calculations
imap_te         = 1; % imap_te=1  ! mapping of the bootstrap coefficients w.r.t Te,
Nzp             = 0; % N from Landremans Gbar/(i-N)=(G+NI)/(i-N)
UseBfore        = 1; % 1: epsilon = (Bmax - Bmin) / (Bmax + Bmin);
                     % 0: epsilon = Rmin/Rc %Rc =(R_+ + R_-)/2 see IIIA of Hager & Chang (2016)
                     % where R_+ and R_- are major radii at the outer and inner-most points of a flux surface         
filename = 'neo_input.nc'; % neo_input.nc filename

if(imap_te==1)
    % QA_Fig1 ~/src/fusion-io/util/_stellar/write_neo_input -bootstrap 0 -m3dc1 ../C1.h5 -1 -te_start 11950 -te_end 1500 -nr 40 -nphi 540 -ntheta 400 -tol 1 -R0 12.27
    % set according to the write_neo_input command used to generate neo_input.nc file
    Te_psi_avgtart = 11950;
    te_end         = 1500;
else
    % Circ1 ~/src/fusion-io/util/_stellar/write_neo_input -bootstrap 0 -m3dc1 ../C1.h5 -1 -psi_start 0.70 -psi_end 0.96 -nr 100 -tol 0.001
    % set psi_start and psi_end according to the write_neo_input command used to generate neo_input.nc file
    psi_start = 0.70;
    psi_end   = 0.96;
end

nfp=1             % number of field periods set as 1, no change needed for this parameter as of now

%%
%========================================================
% Load data from netCDF file
%========================================================
[data,nr,nphi,ntheta] = read_neo_input(filename);
q_of_psi=abs(data.q);

%%
%========================================================
%Defining Constants
%========================================================
psi_count=nr;
lambda_count=100;

istart=1;
jstart=1;

tol=0.001;
Zcharge=1;
Zcharge_eff=1;

if(imap_te==1)
    psi_norm = (data.temax-data.Te)/data.temax;
else
    for s=1:nr
        if(s==1)
	        psi_norm (s,1)= psi_start;
        else 
	        psi_norm (s,1)= (psi_end-psi_start)*(s-1)/(nr-1.) + psi_start;
        end
    end
end


value_psin=zeros(nr,nphi,ntheta);
for m=1:nr
    value_psin(m,:,:)=psi_norm(m,1);
end


%%
%========================================================
% Calculate dTe/dPsi
%========================================================
tnorm=1000.0;
for i=1:nr
    if (i==1)
        dtebydpsit(i,1)=(data.Te(i+1,1)-data.Te(i,1))/(data.psit(i+1,1)-data.psit(i,1))/tnorm;
    elseif(i==nr)
        dtebydpsit(i,1)=(data.Te(i,1)-data.Te(i-1,1))/(data.psit(i,1)-data.psit(i-1,1))/tnorm;
    else
        dtebydpsit(i,1)=(data.Te(i+1,1)-data.Te(i-1,1))/(data.psit(i+1,1)-data.psit(i-1,1))/tnorm;
    end
end

%%
%========================================================
%calculating jacobian offline and defining theta and phi grid
%========================================================
%calculating jacobian
jac = jacobian3(nr, nphi, ntheta, data.R,data.Z,nfp);

for i = 1: ntheta
  theta(1,i) = 2*pi*double(i-1)/double(ntheta);
end 

for i = 1: nphi
  phi(1,i) = 2*pi*double(i-1)/double(nphi)/nfp;
end 
%========================================================

%%
%========================================================
%Finding flux surface avg using FSA function
%========================================================
%Defining arrays

Bmax_fpsi=zeros(psi_count,1);  %Bmax is the max B on a given flux surf
Bmin_fpsi=zeros(psi_count,1);  %Bmax is the min B on a given flux surf
for isurf =1:psi_count
    Bmax_fpsi(isurf,1)=max(max(data.Bmag(isurf,:,:)));
    Bmin_fpsi(isurf,1)=min(min(data.Bmag(isurf,:,:)));
end  



%%
%=======================================================
%Landreman's isomorphism parameters
%=======================================================
[Izp_sum, Izt_sum, GplusiI_fa, Gbar_by_iminusN] = calculate_current_integrals(data, nphi, ntheta, nfp, Nzp, nr, istart, jstart, tol);
 


%%
%========================================================
%Calculate Ftrap
%========================================================
ftrap = calc_ftrap(data, psi_count, lambda_count, Bmax_fpsi, psi_norm, value_psin, nphi, ntheta, istart, jstart, tol);
%%
%=======================================================
% Calculating collision frequency & Redl/Sauter terms 
%=======================================================
[Zcharge_e, Zcharge_i, Zeff] = calculate_Zcharge(data.ne, data.ni, psi_count, Zcharge);

Bmag_inversefa=1./data.Bmagfa;
[nu_e_star, nu_i_star, epsilon] = calculate_collision_frequencies(psi_count, psi_norm, ...
    data.ne, data.ni, data.Te, data.Ti, Zcharge, Zcharge_i, ...
    Bmax_fpsi, Bmin_fpsi, q_of_psi, data.R, value_psin, ...
    nphi, ntheta,  tol, UseBfore, GplusiI_fa, Bmag_inversefa, Nzp);
    
if UseSauter0Redl1==1
    [f_t31_1D, f_t32_ee_1D, f_t32_ei_1D, f_t33_1D, f_t34_1D, alpha_nu_i_1D, ...
        L31_1D, L32_1D, L34_1D, F32ee_1D, F32ei_1D] = calculate_redl_terms(psi_count, Zeff, ...
        ftrap, nu_e_star, nu_i_star, Zcharge_eff);
elseif UseSauter0Redl1==0
    [f_t31_1D, f_t32_ee_1D, f_t32_ei_1D, f_t33_1D, f_t34_1D, alpha_nu_i_1D,...
    L31_1D, L32_1D, L34_1D, F32ee_1D, F32ei_1D] = calculate_sauter_terms(psi_count, Zeff, ...
    ftrap, nu_e_star, nu_i_star);
end

%%
%=======================================================
% Writing Output Profile files
%=======================================================
write_profile_files(imap_te,psi_count, L31_1D, L32_1D, L34_1D, ...
    alpha_nu_i_1D, dtebydpsit, Gbar_by_iminusN, ftrap, ...
    Bmax_fpsi, data.B2fa, data.Bxfa, data.Byfa, data.Bzfa, psi_norm, data.ne, data.ni, data.Te, data.Ti);


%%
% %=======================================================
% % Plots for q, ftrap, collisionality, L31,L32,L34,aplha, temp, density,
% % 1/<B^2>, dte/dpsit
% %=======================================================
% 
% 
% %========================================================
%     % Defaults for plots
% width = 3;        % Width in inches 
% height = 2.25;    % Height in inches 
% alw = 0.75;       % AxesLineWidth 
% fsz = 12;         % Fontsize 
% lw = 2;           % LineWidth 
% msz = 8;          % MarkerSize 
% %=======================================================
% 
%     if imap_te == 1
%         xdata =  data.Te;
%         xlabeldata='$Te$'
%     else
%         xdata =  psi_norm;
%         xlabeldata='$\psi$'
%     end
% 
%     % Plot q 
%     plot_generic(xdata, q_of_psi, xlabeldata, 'q', '$q(\psi)$', 'qvste', width, height, fsz, lw, alw, msz);
% 
%     % Plot ftrap 
%     plot_generic(xdata, ftrap, xlabeldata, 'ftrap', '$f_{trap}$', 'ftrap_trapz', width, height, fsz, lw, alw, msz);
% 
%     % Plot Collisionality
%     figure;
%     plot(xdata, alpha_nu_i_1D, 'LineStyle', '--', 'Color', [0 0.4 0.8], 'LineWidth', lw);
%     title('Collisionality');
%     xlabel(xlabeldata, 'Interpreter', 'latex', 'FontSize', fsz);
%     ylabel('$\alpha$', 'Interpreter', 'latex', 'FontSize', fsz);
%     legend('$\alpha$', 'Interpreter', 'latex', 'LineWidth', 1, 'location', 'northeast');
%     set_plot_properties(fsz, lw, alw, width, height ,msz);
%     print('Coeffs_alphavspsi_1D', '-dpng', '-r300');
% 
%     % Plot Coefficients (L31, L32, L34)
%     figure;
%     plot(xdata, L31_1D, 'LineStyle', '-', 'Color', [0 0 1], 'LineWidth', lw); hold on;
%     plot(xdata, L32_1D, 'LineStyle', ':', 'Color', [1 0 0], 'LineWidth', lw);
%     plot(xdata, L34_1D, 'LineStyle', '--', 'Color', [0 0.4 0.8], 'LineWidth', lw); hold off;
%     xlabel(xlabeldata, 'Interpreter', 'latex', 'FontSize', fsz);
%     legend('$L_{31}$', '$L_{32}$', '$L_{34}$', 'Interpreter', 'latex', 'LineWidth', 1, 'location', 'northwest');
%     set_plot_properties(fsz, lw, alw, width, height, msz);
%     print('Coeffs_eachvspsi_1D', '-dpng', '-r300');
% 
%     % Plot 1/B^2 vs Te
%     plot_generic(data.Te, 1 ./ data.B2fa, 'Te', '1/B^2', '$1/<B^2>$', '1byB2', width, height, fsz, lw, alw, msz);
% 
% 
%     % Plot Temperature vs psi
%     figure;
%     plot(psi_norm, data.Ti, 'LineStyle', '-', 'Color', [0 0 1], 'LineWidth', lw); hold on;
%     plot(psi_norm, data.Te, 'LineStyle', '--', 'Color', [1 0 0], 'LineWidth', lw); hold off;
%     title('Temperature');
%     xlabel('$\psi$', 'Interpreter', 'latex', 'FontSize', fsz);
%     ylabel('$Temperature \: (eV)$', 'Interpreter', 'latex', 'FontSize', fsz);
%     legend('$T_i$', '$T_e$', 'Interpreter', 'latex', 'LineWidth', 1, 'location', 'southwest');
%     set_plot_properties(fsz, lw, alw, width, height, msz);
%     print('Temp', '-dpng', '-r300');
% 
%     % Plot Density vs psi
%     figure;
%     plot(psi_norm, data.ni, 'LineStyle', '-', 'Color', [0 0 1], 'LineWidth', lw); hold on;
%     plot(psi_norm, data.ne, 'LineStyle', '--', 'Color', [1 0 0], 'LineWidth', lw); hold off;
%     title('Density');
%     xlabel('$\psi$', 'Interpreter', 'latex', 'FontSize', fsz);
%     ylabel('$n$', 'Interpreter', 'latex', 'FontSize', fsz);
%     legend('$n_i$', '$n_e$', 'Interpreter', 'latex', 'LineWidth', 1, 'location', 'northwest');
%     set_plot_properties(fsz, lw, alw, width, height, msz);
%     print('Density', '-dpng', '-r300');
% 
% 
%     % Plot dTe/dpsi_t
%     changetoev = 1.6022e-9 * (4 * pi * 1e14) / (1e4^2); % Constants for conversion
%     plot_generic(psi_norm, dtebydpsit * changetoev, '$\psi$', '$dT_e/d\psi_t$', '$dT_e/d\psi_t$', 'dtebydpsit', width, height, fsz, lw, alw, msz);
% 
% 
% 

%%
%=======================================================
% ---------------------- Functions Below ----------------------
%=======================================================



% Generic Plot Function
function plot_generic(x_data, y_data, xlabel_text, ylabel_text, plot_title, file_name, width, height, fsz, lw, alw, msz)
    figure;
    plot(x_data, y_data, 'LineStyle', '-', 'Color', [0 0 1], 'LineWidth', lw);
    title(plot_title, 'Interpreter', 'latex', 'FontSize', fsz);
    xlabel(xlabel_text, 'Interpreter', 'latex', 'FontSize', fsz);
    ylabel(ylabel_text, 'Interpreter', 'latex', 'FontSize', fsz);
    set_plot_properties(fsz, lw, alw, width, height,msz);
    print(file_name, '-dpng', '-r300');
end

% Function to set common plot properties
function set_plot_properties(fsz, lw, alw, width, height,msz)
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
    set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
    set(0, 'defaultLineLineWidth', lw);   % Set the default line width to lw
    set(0, 'defaultLineMarkerSize', msz); % Set the default line marker size to msz
    set(0, 'defaultFigureInvertHardcopy', 'on'); % This is the default anyway
    set(0, 'defaultFigurePaperUnits', 'inches'); % This is the default anyway
    papersize = get(gcf, 'PaperSize');
    left = (papersize(1) - width) / 2;
    bottom = (papersize(2) - height) / 2;
    myfiguresize = [left, bottom, width, height];
    set(gcf, 'PaperPosition', myfiguresize);
    grid on;
end



%=======================================================
% ---------------------- Read Data ----------------------
function [data, nr, nphi, ntheta] = read_neo_input(filename)
    % Reads the given NETCDF file (filename) and extracts relevant variables.


    % Reading variables from the NetCDF file
    data.temax = ncread(filename, "temax");
    data.psi = ncread(filename, "psi0");
    data.q = ncread(filename, "q");
    data.Te = ncread(filename, "Te0");
    data.Ti = ncread(filename, "Ti0");
    data.ne = ncread(filename, "ne0");
    data.ni = ncread(filename, "ni0");
    data.Br = ncread(filename, "B_R");
    data.Bphi = ncread(filename, "B_Phi");
    data.B_Z = ncread(filename, "B_Z");

    data.Bxfa = ncread(filename, "bx_fa");
    data.Byfa = ncread(filename, "by_fa");
    data.Bzfa = ncread(filename, "bz_fa");
    data.Bmagfa = ncread(filename, "bmag_fa");
    data.B2fa = ncread(filename, "b2_fa");

    data.R = ncread(filename, "R");
    data.Z = ncread(filename, "Z");
    data.Phi = ncread(filename, "Phi");
    data.axis_r = ncread(filename, "axis_R");
    data.axis_z = ncread(filename, "axis_Z");

    data.Te3d = ncread(filename, "te_3d");
    data.Ti3d = ncread(filename, "ti_3d");
    data.ne3d = ncread(filename, "ne_3d");
    data.ni3d = ncread(filename, "ni_3d");

    % Compute the magnetic field magnitude
    data.Bmag = sqrt(data.Br.^2 + data.B_Z.^2 + data.Bphi.^2);

    data.psip = ncread(filename, "Psi_p");
    data.psit = ncread(filename, "Psi_t");

    data.jac_m3dc1 = ncread(filename, "Jac");

    % Get the dimensions of the 3D data arrays
    [nr, nphi, ntheta] = size(data.ne3d);
end


%=======================================================
% ---------------------- Landreman's isomorphism parameters ----------------------
function [Izp_sumj, Izt_sumk, GplusiI_fa, Gbar_by_iminusN] = calculate_current_integrals(data, nphi, ntheta, nfp, Nzp, nr, istart, jstart, tol)

% Loop over all radial, poloidal and toroidal grid points
 for i=1:nr
    Izp_sum(i,1)=0;
    Izt_sum(i,1)=0;
    for j=1:nphi
        for k=1:ntheta

            if(j==1)
                dRdj = (data.R(i,j+1,k) - data.R(i,j,k));
	            dZdj = (data.Z(i,j+1,k) - data.Z(i,j,k));
            elseif (j==nphi)
                dRdj = (data.R(i,j,k) - data.R(i,j-1,k));
	            dZdj = (data.Z(i,j,k) - data.Z(i,j-1,k));
            else
                dRdj = (data.R(i,j+1,k) - data.R(i,j-1,k))/2;
	            dZdj = (data.Z(i,j+1,k) - data.Z(i,j-1,k))/2;
            end
    
    
            if(k==1)
                dRdk = (data.R(i,j,k+1) - data.R(i,j,k));
	            dZdk = (data.Z(i,j,k+1) - data.Z(i,j,k));
            elseif (k==ntheta)
                dRdk = (data.R(i,j,k) - data.R(i,j,k-1));
	            dZdk = (data.Z(i,j,k) - data.Z(i,j,k-1));
            else
                dRdk = (data.R(i,j,k+1) - data.R(i,j,k-1))/2;
	            dZdk = (data.Z(i,j,k+1) - data.Z(i,j,k-1))/2;
            end
            dphidj = 2.* pi / nphi/ nfp;

            %B.dl=d.xj=RBphi phij + Br Rj + BzZj
            Izp_here(i,j,k)=data.R(i,j,k)*data.Bphi(i,j,k)*dphidj+data.Br(i,j,k)*dRdj + data.B_Z(i,j,k)*dZdj;
            Izp_sum(i,1)=Izp_sum(i,1)+Izp_here(i,j,k);

            %B.dl=d.xj= Br Rk + BzZk
            Izt_here(i,j,k)=data.Br(i,j,k)*dRdk + data.B_Z(i,j,k)*dZdk;%+Br(i,j,k)*dRdi + B_Z(i,j,k)*dZdi;
            Izt_sum(i,1)=Izt_sum(i,1)+Izt_here(i,j,k);
        end
    end
  
end

for i=1:nr
    Izt_sumk(i,1)=(Izt_sum(i,1))/nphi/(2*pi);
    Izp_sumj(i,1)=(Izp_sum(i,1))/ntheta/(2*pi);
end

% Calculate GplusiI_fa and Gbar_by_iminusN
GplusiI_fa(:,1)=Izp_sumj(:,1)+Izt_sumk(:,1)./data.q(:,1);

Gbar_by_iminusN=(Izp_sumj(:,1)+Nzp*Izt_sumk(:,1))./(1./data.q(:,1)-Nzp);


end


%=======================================================
% ---------------------- trapped fraction ----------------------
function ftrap = calc_ftrap(data, psi_count, lambda_count, Bmax_fpsi, psi_eval, value_psin, nphi, ntheta, istart, jstart, tol)
%========================================================
%========================================================
%Finding f_trap as a function of psi
%F_𝑡𝑟𝑎𝑝=1−𝑓_𝑒=1−3/4 ⟨𝐵^2 ⟩ ∫_0^(〖1/𝐵〗_𝑚𝑎𝑥)▒𝜆𝑑𝜆/⟨√(1−𝜆𝐵)⟩ 
%========================================================
%========================================================
%========================================================
    Bmag=data.Bmag;
    B2_psi_avg=data.B2fa;
    jac_m3dc1=data.jac_m3dc1;
    % Define dlambda and lambda for every flux surface (0-1/Bmax)
    for m =1:psi_count
        dl(m,1)=(1/Bmax_fpsi(m,1)-0)/(lambda_count);  
        lambda(m,1)=0+dl(m,1)/2;
        for n = 2:lambda_count
            lambda(m,n)=lambda(m,1)+dl(m,1)*(n-1);
        end
    end
    
    %changing the limits of lambda from 0 to lambda(lambda_count) to
    % 1/2(b-a)u+1/2(b+a)
    gauss_count_n=5;
    if(gauss_count_n==5)
        abcissae(1,1)=-0.90618;
        abcissae(2,1)=-0.53847;
        abcissae(3,1)=0.00000;
        abcissae(4,1)=0.53847;
        abcissae(5,1)=0.9061;
        
        weight(1,1)=0.23693;
        weight(2,1)=0.47863;
        weight(3,1)=0.56889;
        weight(4,1)=0.47863;
        weight(5,1)=0.23693;
    end
    
    for i=1:gauss_count_n
        u(i,1)=abcissae(i,1);
    end
    
    %defining dlambda and lambda for every flux surface (0-1/Bmax)
    for m =1:psi_count 
        for i=1:gauss_count_n
            lambda_gauss(m,i)= 1/2*(lambda(m,lambda_count)-lambda(m,1))*u(i,1)+1/2*(lambda(m,lambda_count)+lambda(m,1));
        end
    end

    denom_gauss=zeros(psi_count,gauss_count_n);
    counter3=zeros(psi_count,gauss_count_n);
    for m =1:psi_count  
        current_psi=psi_eval(m,1);   
        for n =1:gauss_count_n
           counter3(m,n)=0;
           counter_negativesqrt_gauss(m,n)=0;
           dV=0;
           for i =istart:nphi
              for j =jstart:ntheta
                 if (value_psin(m,i,j)<=current_psi+tol && value_psin(m,i,j)>=current_psi-tol)
                    if(lambda_gauss(m,n)*Bmag(m,i,j)>1 )
                       %if value inside sqrt is negative set it to
                       %0 and increase the counter for that
                       %particular lambda and psi, generally it
                       %would happen on the last value of lamda
                       %for a given psi
                       if(abs(lambda_gauss(m,n)*Bmag(m,i,j)-1)<0.00001)
                          % skip if negative
                       else
                          counter_negativesqrt_gauss(m,n)=counter_negativesqrt_gauss(m,n)+1;
                       end
                       counter3(m,n)=counter3(m,n)+1;
                     else
                       denom_gauss(m,n)=denom_gauss(m,n)+sqrt(1-lambda_gauss(m,n)*Bmag(m,i,j))*jac_m3dc1(m,i,j);
                       counter3(m,n)=counter3(m,n)+1;
                       dV=dV+jac_m3dc1(m,i,j);
                    end  
                 end
              end
           end    
           %denominator is a flux surface average of sqrt(1-lambda B)
           denom_gauss(m,n)=denom_gauss(m,n)/dV;
           f_gauss(m,n)=lambda_gauss(m,n)/denom_gauss(m,n);
        end              
    end
    
    integral_gauss=zeros(psi_count,1);
    for m =1:psi_count
         %finding integral using trapezoid rule
        for n=1:gauss_count_n
            integral_gauss(m,1)=integral_gauss(m,1)+1/2*(lambda(m,lambda_count)-lambda(m,1))*weight(n,1)*(f_gauss(m,n));
        end
    end
    
    
    for m =1:psi_count
        ftrap_gaussian(m,1)=1-3/4*B2_psi_avg(m,1)* integral_gauss(m,1); 
    end

    ftrap=ftrap_gaussian;
end


%=======================================================
% ---------------------- Zcharge ----------------------
function [Zcharge_e, Zcharge_i, Zeff] = calculate_Zcharge(ne_psi_avg, ni_psi_avg, psi_count, Zcharge)
    % Constants
    Zs = 1;           % Charge state for electrons
    Z_alpha = 1;      % Parameter for ion charge calculation
    
    % Initialize total sums
    ne_psi_tot = 0;
    ni_psi_tot = 0;
    
    % Sum over all psi
    for m = 1:psi_count  
        ne_psi_tot = ne_psi_tot + ne_psi_avg(m, 1);
        ni_psi_tot = ni_psi_tot + ni_psi_avg(m, 1);
    end  
    
    % Calculate Z_bar and Zeff
    Z_bar = ne_psi_tot / ni_psi_tot;
    Zeff = (ni_psi_tot) * Zs^2 / ne_psi_tot;
    
    % Calculate Zcharge_i
    Zcharge_e = Zcharge;  % Zcharge for electrons (given as input)
    Zcharge_i = (Z_alpha^2 * Z_bar * Zeff) ^ 0.25;  % Zcharge for ions
    
end


%=======================================================
% ---------------------- Collision frequency ----------------------
function [nu_e_star, nu_i_star, epsilon] = calculate_collision_frequencies(psi_count, psi_eval, ...
    ne_psi_avg, ni_psi_avg, Te_psi_avg, Ti_psi_avg, Zcharge, Zcharge_i, ...
    Bmax_fpsi, Bmin_fpsi, q_of_psi, R, value_psin, ...
    nphi, ntheta,  tol, UseBfore, GplusiI_fa, Bmag_inversefa, Nzp)
    
    % Initialize output arrays
    nu_e_star = zeros(psi_count, 1);
    nu_i_star = zeros(psi_count, 1);
    ln_lambda_e = zeros(psi_count, 1);
    ln_lambda_i = zeros(psi_count, 1);
    epsilon = zeros(psi_count, 1);
    qR = zeros(psi_count, 1);
    
    % Loop over all flux surfaces (psi_count)
    for m = 1:psi_count
        current_psi = psi_eval(m, 1);
        
        % Calculate epsilon as a function of psi (based on whether UseBfore is set)
        if UseBfore == 1
            % Using Landreman et al (2022) epsilon calculation
            epsilon(m, 1) = (Bmax_fpsi(m, 1) - Bmin_fpsi(m, 1)) / (Bmax_fpsi(m, 1) + Bmin_fpsi(m, 1));
            
            % Calculate qR using a specified method 
            %qR= (G+iI)/(i-N) <1/B>
            qR_orig = GplusiI_fa(m, 1) * Bmag_inversefa(m, 1) / (1 / q_of_psi(m, 1) - Nzp);
            %qR_usingjac = jacB2_fa(m, 1) * Bmag_inversefa(m, 1) / (1 / q_of_psi(m, 1) - Nzp);
            %qR_usingRfa = Rfa(m, 1) * q_of_psi(m, 1);
            qR(m, 1) = qR_orig;  % Using the original qR method (or any other if necessary)
            
            % Calculate ln_lambda_e and nu_e_star for electrons
            ln_lambda_e(m, 1) = 31.3 - log(sqrt(ne_psi_avg(m, 1)) / Te_psi_avg(m, 1));
            temp = 6.921 * 10^(-18) * Zcharge * ne_psi_avg(m, 1) * ln_lambda_e(m, 1) / Te_psi_avg(m, 1)^2 / epsilon(m, 1)^(3 / 2) * qR(m, 1);
            nu_e_star(m, 1) = abs(temp);
            
            % Calculate ln_lambda_i and nu_i_star for ions
            Zcharge = Zcharge_i;  % Use ion charge state
            ln_lambda_i(m, 1) = 30 - log(Zcharge^3 * sqrt(ni_psi_avg(m, 1)) / Ti_psi_avg(m, 1)^(3 / 2));
            temp = 4.9 * 10^(-18) * ni_psi_avg(m, 1) * Zcharge^4 * ln_lambda_i(m, 1) / Ti_psi_avg(m, 1)^2 / epsilon(m, 1)^(3 / 2) * qR(m, 1);
            nu_i_star(m, 1) = abs(temp);
        else
            % Alternative epsilon calculation when UseBfore is not used
            %mean minor radius Rmin=(R_+ - R_-)/2 see IIIA of Hager & Chang (2016)
            %Geometrical Center Rc =(R_+ + R_-)/2 see IIIA of Hager & Chang (2016)
            % where R_+ and R_- are major radii at the outer and inner-most points of a flux surface 
            % epsilon=Rmin/Rc
            rad = zeros(nphi, ntheta);
            counter = 0;
            for i = 1:nphi
                for j = 1:ntheta
                    if (value_psin(m, i, j) <= current_psi + tol / 6 && value_psin(m, i, j) >= current_psi - tol / 6)
                        rad(i, j) = R(m, i, j);
                        counter = counter + 1;
                    end
                end
            end
            
            % Calculate R_plus, R_minus, and epsilon
            R_plus = max(max(rad));
            rad(rad == 0) = inf;
            R_minus = min(min(rad));
            
            Rmin = (R_plus - R_minus) / 2;
            Rc = (R_plus + R_minus) / 2;
            epsilon(m, 1) = Rmin / Rc;
            
            % Calculate qR using the radius of curvature method
            qR(m, 1) = q_of_psi(m, 1) * Rc;
            
            % Calculate ln_lambda_e and nu_e_star for electrons
            ln_lambda_e(m, 1) = 31.3 - log(sqrt(ne_psi_avg(m, 1)) / Te_psi_avg(m, 1));
            temp = 6.921 * 10^(-18) * ne_psi_avg(m, 1) * ln_lambda_e(m, 1) / Te_psi_avg(m, 1)^2 / epsilon(m, 1)^(3 / 2) * q_of_psi(m, 1) * Rc;
            nu_e_star(m, 1) = temp;
            
            % Calculate ln_lambda_i and nu_i_star for ions
            Zcharge = Zcharge_i;  % Use ion charge state
            ln_lambda_i(m, 1) = 30 - log(Zcharge^3 * sqrt(ni_psi_avg(m, 1)) / Ti_psi_avg(m, 1)^(3 / 2));
            temp = 4.9 * 10^(-18) * ni_psi_avg(m, 1) * Zcharge^4 * ln_lambda_i(m, 1) / Ti_psi_avg(m, 1)^2 / epsilon(m, 1)^(3 / 2) * q_of_psi(m, 1) * Rc;
            nu_i_star(m, 1) = temp;
        end
    end
end

%=======================================================
% ---------------------- Redl Paramters ----------------------
function [f_t31_1D, f_t32_ee_1D, f_t32_ei_1D, f_t33_1D, f_t34_1D, alpha_nu_i_1D, ...
    L31_1D, L32_1D, L34_1D, F32ee_1D, F32ei_1D] = calculate_redl_terms(psi_count, Zeff, ...
    ftrap, nu_e_star, nu_i_star, Zcharge_eff)

    % Initialize output arrays
    f_t31_1D = zeros(psi_count, 1);
    f_t32_ee_1D = zeros(psi_count, 1);
    f_t32_ei_1D = zeros(psi_count, 1);
    f_t33_1D = zeros(psi_count, 1);
    f_t34_1D = zeros(psi_count, 1);
    alpha_nu_i_1D = zeros(psi_count, 1);

    L31_1D = zeros(psi_count, 1);
    L32_1D = zeros(psi_count, 1);
    L34_1D = zeros(psi_count, 1);
    F32ee_1D = zeros(psi_count, 1);
    F32ei_1D = zeros(psi_count, 1);

    % Loop through the psi_count to calculate terms for each value of psi
    for m = 1:psi_count  
        
        % Calculate f_t31_1D
        temp = 1 + ...
               0.67 * (1 - 0.7 * ftrap(m, 1)) * sqrt(nu_e_star(m, 1)) / (0.56 + 0.44 * Zcharge_eff) + ...
               (0.52 + 0.086 * sqrt(nu_e_star(m, 1))) * (1 + 0.87 * ftrap(m, 1)) * nu_e_star(m, 1) / ...
               (1 + 1.13 * (Zcharge_eff - 1)^0.5);
        f_t31_1D(m, 1) = ftrap(m, 1) / temp;

        % Calculate f_t32_ee_1D
        temp = 1 + ...
               0.23 * (1 - 0.96 * ftrap(m, 1)) * sqrt(nu_e_star(m, 1)) / Zcharge_eff^0.5 + ...
               0.13 * (1 - 0.38 * ftrap(m, 1)) * nu_e_star(m, 1) / Zcharge_eff^2 * ...
               (sqrt(1 + 2 * (Zcharge_eff - 1)^0.5) + ftrap(m, 1)^2 * sqrt((0.075 + 0.25 * (Zcharge_eff - 1)^2)) * nu_e_star(m, 1));
        f_t32_ee_1D(m, 1) = ftrap(m, 1) / temp;

        % Calculate f_t32_ei_1D
        temp = 1 + ...
               0.87 * (1 + 0.39 * ftrap(m, 1)) * sqrt(nu_e_star(m, 1)) / (1 + 2.95 * (Zcharge_eff - 1)^2) + ...
               1.53 * (1 - 0.37 * ftrap(m, 1)) * nu_e_star(m, 1) * (2 + 0.375 * (Zcharge_eff - 1));
        f_t32_ei_1D(m, 1) = ftrap(m, 1) / temp;

        % Calculate f_t33_1D
        temp = 1 + ...
               0.25 * (1 - 0.7 * ftrap(m, 1)) * sqrt(nu_e_star(m, 1)) * (1 + 0.45 * (Zcharge_eff - 1)^0.5) + ...
               (0.61 * (1 - 0.41 * ftrap(m, 1)) * nu_e_star(m, 1)) / Zcharge_eff * 0.5;
        f_t33_1D(m, 1) = ftrap(m, 1) / temp;

        % Calculate alpha_nu_i_1D
        alpha_0 = -(0.62 + 0.055 * (Zcharge_eff - 1)) / (0.53 + 0.17 * (Zcharge_eff - 1)) * ...
                  (1 - ftrap(m, 1)) / (1 - (0.31 - 0.065 * (Zcharge_eff - 1)) * ftrap(m, 1) - 0.25 * ftrap(m, 1)^2);
        alpha_nu_i_1D(m, 1) = ((alpha_0 + 0.7 * Zcharge_eff * ftrap(m, 1)^0.5 * sqrt(nu_e_star(m, 1))) / ...
                               (1 + 0.18 * sqrt(nu_i_star(m, 1))) - 0.002 * nu_i_star(m, 1)^2 * ftrap(m, 1)^6) / ...
                              (1 + 0.004 * nu_i_star(m, 1)^2 * ftrap(m, 1)^6);

        % Calculate L31_1D
        L31_1D(m, 1) = (1 + 0.15 / (Zcharge_eff^1.2 - 0.71)) * f_t31_1D(m, 1) - ...
                        0.22 / (Zcharge_eff^1.2 - 0.71) * f_t31_1D(m, 1)^2 + ...
                        0.01 / (Zcharge_eff^1.2 - 0.71) * f_t31_1D(m, 1)^3 + ...
                        0.06 / (Zcharge_eff^1.2 - 0.71) * f_t31_1D(m, 1)^4;

        % Calculate L34_1D
        L34_1D(m, 1) = L31_1D(m, 1);

        % Calculate F32ee_1D
        temp = f_t32_ee_1D(m, 1);
        F32ee_1D(m, 1) = (0.1 + 0.6 * Zcharge_eff) / (Zcharge_eff * (0.77 + 0.63 * (1 + (Zcharge_eff - 1)^1.1))) * ...
                          (temp - temp^4) + ...
                          (0.7) / (1 + 0.2 * Zcharge_eff) * (temp^2 - temp^4 - 1.2 * (temp^3 - temp^4)) + ...
                          (1.3) / (1 + 0.5 * Zcharge_eff) * temp^4;

        % Calculate F32ei_1D
        temp = f_t32_ei_1D(m, 1);
        F32ei_1D(m, 1) = -(0.4 + 1.93 * Zcharge_eff) / (Zcharge_eff * (0.8 + 0.6 * Zcharge_eff)) * ...
                          (temp - temp^4) + ...
                          5.5 / (1.5 + 2 * Zcharge_eff) * (temp^2 - temp^4 - 0.8 * (temp^3 - temp^4)) - ...
                          (1.3) / (1 + 0.5 * Zcharge_eff) * temp^4;

        % Calculate L32_1D
        L32_1D(m, 1) = F32ee_1D(m, 1) + F32ei_1D(m, 1);
    end
end


%=======================================================
% ---------------------- Sauter Paramters ----------------------


function [f_t31_1D, f_t32_ee_1D, f_t32_ei_1D, f_t33_1D, f_t34_1D, alpha_nu_i_1D,...
    L31_1D, L32_1D, L34_1D, F32ee_1D, F32ei_1D] = calculate_sauter_terms(psi_count, Zeff, ...
    ftrap, nu_e_star, nu_i_star)
    % Initialize output arrays
    f_t31_1D = zeros(psi_count, 1);
    f_t32_ee_1D = zeros(psi_count, 1);
    f_t32_ei_1D = zeros(psi_count, 1);
    f_t33_1D = zeros(psi_count, 1);
    f_t34_1D = zeros(psi_count, 1);
    alpha_nu_i_1D = zeros(psi_count, 1);

    L31_1D = zeros(psi_count, 1);
    L32_1D = zeros(psi_count, 1);
    L34_1D = zeros(psi_count, 1);
    F32ee_1D = zeros(psi_count, 1);
    F32ei_1D = zeros(psi_count, 1);

    % Loop over psi_count and perform the calculations
    for m = 1:psi_count
        % Calculate f_t31_1D
        Zcharge = Zeff;
        temp = 1 + (1 - 0.1 * ftrap(m, 1)) * sqrt(nu_e_star(m, 1)) + 0.5 * (1 - ftrap(m, 1)) * nu_e_star(m, 1) / Zcharge;
        f_t31_1D(m, 1) = ftrap(m, 1) / temp;

        % Calculate f_t32_ee_1D
        Zcharge = Zeff;
        temp = 1 + 0.26 * (1 - ftrap(m, 1)) * sqrt(nu_e_star(m, 1)) + 0.18 * (1 - 0.37 * ftrap(m, 1)) * nu_e_star(m, 1) / Zcharge^0.5;
        f_t32_ee_1D(m, 1) = ftrap(m, 1) / temp;

        % Calculate f_t32_ei_1D
        Zcharge = Zeff;
        temp = 1 + (1 + 0.6 * ftrap(m, 1)) * sqrt(nu_e_star(m, 1)) + 0.85 * (1 - 0.37 * ftrap(m, 1)) * nu_e_star(m, 1) * (1 + Zcharge);
        f_t32_ei_1D(m, 1) = ftrap(m, 1) / temp;

        % Calculate f_t34_1D
        Zcharge = Zeff;
        temp = 1 + (1 - 0.1 * ftrap(m, 1)) * sqrt(nu_e_star(m, 1)) + 0.5 * (1 - 0.5 * ftrap(m, 1)) * nu_e_star(m, 1) / Zcharge;
        f_t34_1D(m, 1) = ftrap(m, 1) / temp;

        % Calculate alpha_nu_i_1D
        alpha_0 = -1.17 * (1 - ftrap(m, 1)) / (1 - 0.22 * ftrap(m, 1) - 0.19 * ftrap(m, 1)^2);
        alpha_nu_i_1D(m, 1) = ((alpha_0 + 0.25 * (1 - ftrap(m, 1)^2) * sqrt(nu_i_star(m, 1))) / (1 + 0.5 * sqrt(nu_i_star(m, 1))) + 0.315 * nu_i_star(m, 1)^2 * ftrap(m, 1)^6) * 1 / (1 + 0.15 * nu_i_star(m, 1)^2 * ftrap(m, 1)^6);

        % Calculate L31_1D
        Zcharge = Zeff;
        L31_1D(m, 1) = (1 + 1.4 / (Zcharge + 1)) * f_t31_1D(m, 1) ...
            - 1.9 / (Zcharge + 1) * f_t31_1D(m, 1)^2 ...
            + 0.3 / (Zcharge + 1) * f_t31_1D(m, 1)^3 ...
            + 0.2 / (Zcharge + 1) * f_t31_1D(m, 1)^4;

        % Calculate L34_1D
        Zcharge = Zeff;
        L34_1D(m, 1) = (1 + 1.4 / (Zcharge + 1)) * f_t34_1D(m, 1) ...
            - 1.9 / (Zcharge + 1) * f_t34_1D(m, 1)^2 ...
            + 0.3 / (Zcharge + 1) * f_t34_1D(m, 1)^3 ...
            + 0.2 / (Zcharge + 1) * f_t34_1D(m, 1)^4;

        % Calculate F32ee_1D
        Zcharge = Zeff;
        temp = f_t32_ee_1D(m, 1);
        F32ee_1D(m, 1) = (0.05 + 0.62 * Zcharge) / (Zcharge * (1 + 0.44 * Zcharge)) * (temp - temp^4) ...
            + (1) / (1 + 0.22 * Zcharge) * (temp^2 - temp^4 - 1.2 * (temp^3 - temp^4)) ...
            + (1.2) / (1 + 0.5 * Zcharge) * temp^4;

        % Calculate F32ei_1D
        Zcharge = Zeff;
        temp = f_t32_ei_1D(m, 1);
        F32ei_1D(m, 1) = -(0.56 + 1.93 * Zcharge) / (Zcharge * (1 + 0.44 * Zcharge)) * (temp - temp^4) ...
            + 4.95 / (1 + 2.48 * Zcharge) * (temp^2 - temp^4 - 0.55 * (temp^3 - temp^4)) ...
            - (1.2) / (1 + 0.5 * Zcharge) * temp^4;

        % Calculate L32_1D
        L32_1D(m, 1) = F32ee_1D(m, 1) + F32ei_1D(m, 1);
    end
end


%=======================================================
% ---------------------- Writing Profiles ----------------------
function write_profile_files(imap_te, psi_count, L31_1D, L32_1D, L34_1D, ...
    alpha_nu_i_1D, dtebydpsit, Gbar_by_iminusN, ftrap, ...
    Bmax_fpsi, B2fa, Bxfa, Byfa, Bzfa, psi_eval, ne_psi_avg, ni_psi_avg, Te_psi_avg, Ti_psi_avg)

    % Write the first file: ProfileJBSCoeff_Te_L31_32_34_alpha_B2_dtedpsit_G
    
    if imap_te == 1
        fileID = fopen('ProfileJBSCoeff_Te_L31_32_34_alpha_B2_dtedpsit_G', 'w');
        if fileID == -1
            error('Could not open ProfileJBSCoeff_Te_L31_32_34_alpha_B2_dtedpsit_G for writing');
        end
        fprintf(fileID, '%9s %9s %9s %9s %9s %9s %9s %9s\n', ...
                'Te', 'L31', 'L32', 'L34', 'alpha', '1/<B^2>', 'dtebydpsit*2*pi', 'Gbar/(i-N)');
        
        for i = 2:psi_count   
            fprintf(fileID, '%6.8f %6.8f %6.8f %6.8f %6.8f %6.8f %6.8f %6.8f\n', ...
                    Te_psi_avg(i, 1), L31_1D(i, 1), L32_1D(i, 1), L34_1D(i, 1), ...
                    alpha_nu_i_1D(i, 1), 1/B2fa(i, 1), dtebydpsit(i, 1) * 2 * pi, Gbar_by_iminusN(i, 1));
        end
        fclose(fileID);
    else
        fileID = fopen('ProfileJBSCoeff_Psi_L31_32_34_alpha_B2_dtedpsit_G', 'w');
        if fileID == -1
            error('Could not open ProfileJBSCoeff_Te_L31_32_34_alpha_B2_dtedpsit_G for writing');
        end
        fprintf(fileID, '%9s %9s %9s %9s %9s %9s %9s %9s\n', ...
                'Psi', 'L31', 'L32', 'L34', 'alpha', '1/<B^2>', 'dtebydpsit*2*pi', 'Gbar/(i-N)');
        
        for i = 2:psi_count   
            fprintf(fileID, '%6.8f %6.8f %6.8f %6.8f %6.8f %6.8f %6.8f %6.8f\n', ...
                    psi_eval(i, 1), L31_1D(i, 1), L32_1D(i, 1), L34_1D(i, 1), ...
                    alpha_nu_i_1D(i, 1), 1/B2fa(i, 1), dtebydpsit(i, 1) * 2 * pi, Gbar_by_iminusN(i, 1));
        end
        fclose(fileID);
    end

    % Write the second file: Ftrap_Bmax_writeneo
    fileID = fopen('Ftrap_Bmax_writeneo', 'w');
    if fileID == -1
        error('Could not open Ftrap_Bmax_writeneo for writing');
    end
    fprintf(fileID, '%9s %9s %9s %9s %9s %9s %9s\n', ...
            'psi', 'ftrap', 'Bmax', 'B^2', 'Br', 'Bphi', 'Bz');
    
    for i = 1:psi_count  
        fprintf(fileID, '%6.8f %6.8f %6.8f %6.8f %6.8f %6.8f %6.8f\n', ...
                psi_eval(i, 1), ftrap(i, 1), Bmax_fpsi(i, 1), B2fa(i, 1), ...
                Bxfa(i, 1), Byfa(i, 1), Bzfa(i, 1));
    end
    fclose(fileID);

    % Write the third file: NeNiTeTi_writeneo
    fileID = fopen('NeNiTeTi_writeneo', 'w');
    if fileID == -1
        error('Could not open NeNiTeTi_writeneo for writing');
    end
    fprintf(fileID, '%9s %9s %9s %9s %9s\n', 'psi', 'ne', 'ni', 'Te', 'Ti');
    
    for i = 1:psi_count  
        fprintf(fileID, '%6.8f %6.8f %6.8f %6.8f %6.8f\n', ...
                psi_eval(i, 1), ne_psi_avg(i, 1), ni_psi_avg(i, 1), ...
                Te_psi_avg(i, 1), Ti_psi_avg(i, 1));
    end
    fclose(fileID);

end





