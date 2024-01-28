%Explore mock x-band data

addpath(genpath('C:\Users\josep\Desktop\PHYC40900_Project TP'));
%addpath(genpath('M:\generalMatlab'));
%addpath(genpath('M:\CurrentInversion\NTNU_SIO_DoppVis\Code'));
%addpath(genpath('/Users/susannestole-hentschel/Documents/projects/Ben/Code'));

%%

dataDir = 'C:\Users\josep\Desktop\PHYC40900_Project TP\Data_NEW';%Path to the cubes
saveRes = 0;%Whether or not to write the results to file. (1- saves it, 0 means it doesnt save)
listing = dir(fullfile(dataDir,'shearing*'));
%% 

fileName ='shearing_curr_res_7.5_dt_1.0_T_1200_U_1.0exp(0.5z)+0.05_psi_0_smax_30_0_surf3d.hdf5';%Filename of cube
etaDataSetName = '/eta';%Dataset to use for cube, either the surface elevation or one of the radar image fields.

info = h5info(fullfile(dataDir,fileName));
dataSetNames = {info.Datasets(:).Name};

if numel(info.Groups)>0
    groupNames = {info.Groups(:).Name};
else
    groupNames = {''};
end

i1 = strfind(fileName,'psi');
i2 = strfind(fileName,'smax');
psi = str2double(fileName(i1+4:i2-2));%Current angle relative to x-axis

eta = h5read(fullfile(dataDir,fileName),etaDataSetName);
x = h5read(fullfile(dataDir,fileName),'/x');
y = h5read(fullfile(dataDir,fileName),'/y');
t = h5read(fullfile(dataDir,fileName),'/t');
z = h5read(fullfile(dataDir,fileName),'/z');
U = h5read(fullfile(dataDir,fileName),'/U');
kT = h5read(fullfile(dataDir,fileName),'/k');
Uk = h5read(fullfile(dataDir,fileName),'/Uk');

dx = mean(diff(x));
dy = mean(diff(y));
dt = mean(diff(t));

savePlots = 0;%Whether or not to save SNR maps and cylindrical cross section plots for each wavenumber.

clear STAT_SA;
verbose = 0;%Whether or not to show the surface images

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make the IMG_SEQ structure. - step 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear IMG_SEQ;

for i = 1:size(eta,3)
    if verbose
        figure(1);
        imagesc(squeeze(eta(:,:,i)));
        colorbar;
        axis image;
        colormap gray;
        %caxis([0,150]);
    end
    %xlim([350,450]);ylim([300,500]);
drawnow;
IMG_SEQ.IMG(:,:,i) = DJIFrame_to_IMG_SEQ_format(squeeze(eta(:,:,i)));
end

[gridY,gridX] = meshgrid(1:size(IMG_SEQ.IMG,2),1:size(IMG_SEQ.IMG,1));
gridX = (gridX-mean(gridX(:)))*dx;
gridY = (gridY-mean(gridY(:)))*dy;

IMG_SEQ.gridX = gridX;%X-Y interchanging due to CopterCurrents convention
IMG_SEQ.gridY = gridY;%X-Y interchanging due to CopterCurrents convention
IMG_SEQ.dx = dy;
IMG_SEQ.dy = dx;
IMG_SEQ.dt = dt;

%The rest of the fields should be unnecessary but populated with dummy variables
%to avoid causing errors.
IMG_SEQ.altitude = 100;
IMG_SEQ.pitch = -90;
IMG_SEQ.roll = 0;
IMG_SEQ.heading = 0;
IMG_SEQ.mediainfo_string = '';

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract the wavenumber-dependent Doppler shift velocities - step 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxFrequency = 4.0;%Approximate max  frequency to consider (converted to max wavenumber below) [rad/sec]
verbose = 1;%Whether or not to show diagnostic figures of SNR for each wavenumber.

dk = 2*pi/(IMG_SEQ.dx*min(size(IMG_SEQ.IMG,1),size(IMG_SEQ.IMG,2)));% wavenumber resolution of spectrum in each spatial window (not strictly true if dx ~= dy, but value only needs to be approximate in practice)
kW = 1*dk;%Half width of wavenumbers bins [rad/m]

%Wavenumber values for Doppler shift extraction.
%wavenumbers = dk*10:kW:maxFrequency^2/9.81;
%wavenumbers = 6*dk:kW:60;
%wavenumbers = 6*kW:4*kW:2.5;
%wavenumbers = 3*kW:kW:0.35;
wavenumbers = 3*dk:dk:0.35;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define parameters specific to the k-dependent Doppler shift extraction and
%PEDM

frequencyLimits = [0.01,20.0];%frequency limits for masking the spectrum [min max], rad/sec

% Ux current limits [m/s]
Ux_limits = 1.0*[-1.0 1.0];

% Uy current limts [m/s]
Uy_limits = 1.0*[-1.0 1.0];

% Current step size [m/s]
U_res = 0.1;

% (optional): whether to include 2nd harmonic of the spectrum in the fit (false by default)
include2ndHarmonic = 0;

% (optional): whether to do the fit in log space (false by default)
logFlag = [];

% (optional) omegaWidthFun: function handle as a function of wavenumber i.e.
%@(k) f(k)...., specifying frequency width of the weighting function in
%frequency-angle space (constant wavenumber). Width is half-width 1/e^2
%point of a Gaussian function.
%omegaWidthFun = @(k) 0.4 + 0.1*k;
omegaWidthFun = @(k) 0.05 + 0.0*k;


%The following OPTIONAL parameters involve post-processing of the Doppler shifts:
%SNR_filter: whether to use a signal-to-noise filter (false by default)
SNR_filter = 0;

%SNR_threshold: threshold signal-to-noise value for above filter (set to 2.0 by default)
SNR_threshold = sqrt(1);

%Peak_filter: whether to use a multiple peaks filter (false by default)
Peak_filter = 0;

%Peak_threshold: peak threshold of maximum value (0.5 by default)
Peak_threshold = 0.5;

%Outlier_filter: whether to use an outlier filter (quartile-based) (false by default)
Outlier_filter = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define the fluid & physical properties: depth, gravitational acceleration,
%surface tension coefficient
Properties = struct('h',Inf,'g',9.81,'T',0.072 / 1000);
[Uym, Uxm] = meshgrid(min(Uy_limits):U_res:max(Uy_limits),min(Ux_limits):U_res:max(Ux_limits));

%Initialize Doppler shift variables
Ux = zeros(1,numel(wavenumbers))*NaN;
Uy = zeros(1,numel(wavenumbers))*NaN;
UxC = zeros(1,numel(wavenumbers))*NaN;
UyC = zeros(1,numel(wavenumbers))*NaN;
UxLS = zeros(1,numel(wavenumbers))*NaN;
UyLS = zeros(1,numel(wavenumbers))*NaN;

SNR_max = zeros(1,numel(wavenumbers))*NaN;

clear DSV;
clear DSVPC;
clear DSVLS;

for jj = 1:numel(wavenumbers)
%%
wavenumberLimits = wavenumbers(jj) + kW*[-1,1];%wavenumber limits for masking the spectrum [min max], rad/sec

STCFIT = struct('Name',{fileName});
STCFIT = generate_STCFIT_for_NSPP(STCFIT,wavenumbers(jj),include2ndHarmonic,logFlag,...
    omegaWidthFun,SNR_filter,SNR_threshold,Peak_filter,Peak_threshold,Outlier_filter);


fit_param = STCFIT.NSPP_fit_param;
fit_param.kWidth = 4*kW;

fit_param.Ux_2D = Uxm;
fit_param.Uy_2D = Uym;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now extract the Doppler shifts 
% 
% tic;
%STCFIT = run_current_fit_depth_profile(IMG_SEQ,STCFIT,windowList);
% toc;

% get spectrum
%Read in the power spectrum.

Spectrum = retrieve_power_spectrum(...
    IMG_SEQ.IMG,....
    IMG_SEQ.dx,...
    IMG_SEQ.dy,...
    IMG_SEQ.dt,wavenumberLimits,frequencyLimits);

if fit_param.include2ndHarmonic
    Spectrum2 = retrieve_power_spectrum(...
        IMG_SEQ.IMG,....
        IMG_SEQ.dx,...
        IMG_SEQ.dy,...
        IMG_SEQ.dt,2*wavenumberLimits,frequencyLimits);
    
    Spectrum.power_Spectrum2 = Spectrum2.power_Spectrum;
    Spectrum.Kx_3D2 = Spectrum2.Kx_3D;
    Spectrum.Ky_3D2 = Spectrum2.Ky_3D;
    Spectrum.W_3D2 = Spectrum2.W_3D;
end

%%
%Get Doppler shift velocities
%tic;
out_DS = get_doppler_shift_velocities_nsp(Spectrum,fit_param,Properties,1);

if savePlots
SNR_fileName = fullfile(SNR_plotDir,...
    sprintf('%04iDeg_%04iRadpkm_SNR',round(latitude*1e4),round(wavenumbers(jj)*1e3)));
print(SNR_fileName,'-dpng');
end
%toc;
DSV(jj) = out_DS;

Ux(jj) = -out_DS.Ux_filt;
Uy(jj) = -out_DS.Uy_filt;
SNR_max(jj) = out_DS.SNR_max;


dtheta = 2*dk / wavenumbers(jj);
Scyl = cylinder_cross_section(Spectrum,dtheta,4);
thetaVals = Scyl.thetaM(1,:);


figure(5);imagesc(Scyl.thetaM(1,:)/pi,Scyl.omegaM(:,1),Scyl.P_k);colorbar;shading flat;set(gca,'YDir','normal');
line(thetaVals/pi,sqrt(1)*sqrt(9.81*wavenumbers(jj)) + 0*thetaVals,'LineStyle',':','Color','r','LineWidth',1.5);
line(thetaVals/pi,sqrt(1)*sqrt(9.81*wavenumbers(jj))+wavenumbers(jj)*(cos(thetaVals)*Ux(jj) + sin(thetaVals)*Uy(jj)),'LineStyle','--','Color','r','LineWidth',1.5);
line(thetaVals/pi,-sqrt(1)*sqrt(9.81*wavenumbers(jj))+wavenumbers(jj)*(cos(thetaVals)*Ux(jj) + sin(thetaVals)*Uy(jj)),'LineStyle','-','Color','r','LineWidth',1.5);

% hold on;plot(thetaVals,sqrt(9.81*wavenumbers(jj))-wavenumbers(jj)*(cos(thetaVals)*Ux(jj) + sin(thetaVals)*Uy(jj)),'--w',');
% hold off;
xlabel('\theta/\pi');ylabel('Frequency [rad/s]');
title(sprintf('$k$: %.2f Rad/m \n$\\Delta\\omega/k$: %.2f m/s',wavenumbers(jj),mean(diff(Scyl.omegaM(:,1))) / wavenumbers(jj)),'FontWeight','normal','Interpreter','latex');
%colormap(flipud(othercolor('PuBu9')));
drawnow;
%%
if savePlots
Spec_fileName = fullfile(polarSpecDir,...
    sprintf('%04iDeg_%04iRadpkm_polarSpec',round(latitude*1e4),round(wavenumbers(jj)*1e3)));
print(Spec_fileName,'-dpng');
end

figure(10);
%yyaxis left;
plot(...
    wavenumbers,Ux,'o',...
    wavenumbers,Uy,'xk',...
    kT,Uk*cosd(psi),'--k',...
    kT,-Uk*sind(psi),'--k',...
    'MarkerSize',6);drawnow;
xlabel('Wavenumber [Rad/m]');ylabel('[m/s]');
ylim(1*[-1,1]);
% yyaxis right;
% plot(wavenumbers,SNR_max,'--');
% ylabel('SNR');
hl = legend('$\tilde{c}_x$','$\tilde{c}_y$','Location','sw');
%hl = legend('$\tilde{c}_x$','$\tilde{c}_y$','$\tilde{c}_x$','$\tilde{c}_y$','SNR','Location','e');

set(hl,'Interpreter','latex');
set(gca,'FontSize',14);
drawnow;





end



UxDSV = -[DSV(:).Ux];
UyDSV = -[DSV(:).Uy];

%Write results to file.
resName = fullfile(dataDir,fileName);
N_dsv = sum(cellfun(@numel, strfind(groupNames,'DSV')));
dsvName = sprintf('/DSV_%02i',N_dsv);
locName = [dsvName,etaDataSetName];

if saveRes
dsName = [locName,'/k'];
dataVar = wavenumbers;
h5create(resName,dsName,size(dataVar));
h5write(resName,dsName,dataVar);

dsName = [locName,'/UxDSV'];
dataVar = UxDSV;
h5create(resName,dsName,size(dataVar));
h5write(resName,dsName,dataVar);

dsName = [locName,'/UyDSV'];
dataVar = UyDSV;
h5create(resName,dsName,size(dataVar));
h5write(resName,dsName,dataVar);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform the PEDM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Run inversion methods

inds = 1:numel(wavenumbers);
out = find_current_depth_profile(wavenumbers(inds),UxDSV(inds),UyDSV(inds),Inf);

figure(1);p = plot(...
    out.PEDM.U1_fun(z),z,'--',...
    out.PEDM.U2_fun(z),z,'-.',...
    UxDSV(inds),-1./(2*wavenumbers(inds)),'s',...
    U*cosd(psi),z,'k',...
    -U*sind(psi),z,'-.k','LineWidth',2);

%ylim([min(z),0]);
%ylim([-1./(2*kVals(inds(1))),0]);
ylim([-1./(2*wavenumbers(1)),0]);
xlim(1.0*[-1,1.0]);
xlabel('[m/s]');ylabel('Depth [m]');

p(3).MarkerFaceColor = 0.6*[1,1,1];
p(3).MarkerEdgeColor = 'k';
p(3).LineWidth = 1.0;

legend('PEDM_x','PEDM_y','EDM_x','U_x(z)','U_y(z)','interpreter','none','location','se');
set(gca,'fontsize',12);
drawnow;

UxE = interp1(kT,Uk(1,:)*cosd(psi),wavenumbers);
UyE = interp1(kT,-Uk(1,:)*sind(psi),wavenumbers);

outExact = find_current_depth_profile(wavenumbers(inds),UxE(inds),UyE(inds),Inf);

if saveRes
    resName = fullfile(dataDir,fileName);
    
    dsName = [locName,'/z'];
    dataVar = z;
    h5create(resName,dsName,size(dataVar));
    h5write(resName,dsName,dataVar);
    
    dsName = [locName,'/UxPEDM'];
    dataVar = out.PEDM.U1_fun(z);
    h5create(resName,dsName,size(dataVar));
    h5write(resName,dsName,dataVar);
        
    dsName = [locName,'/UyPEDM'];
    dataVar = out.PEDM.U2_fun(z);
    h5create(resName,dsName,size(dataVar));
    h5write(resName,dsName,dataVar);
    
    dsName = [locName,'/UxPEDM_exact'];
    dataVar = outExact.PEDM.U1_fun(z);
    h5create(resName,dsName,size(dataVar));
    h5write(resName,dsName,dataVar);
        
    dsName = [locName,'/UyPEDM_exact'];
    dataVar = outExact.PEDM.U2_fun(z);
    h5create(resName,dsName,size(dataVar));
    h5write(resName,dsName,dataVar);
    
end


%%% Repeat not using the first two lowest wavenumbers
inds = 3:numel(wavenumbers);
out = find_current_depth_profile(wavenumbers(inds),UxDSV(inds),UyDSV(inds),Inf);

figure(2);p = plot(...
    out.PEDM.U1_fun(z),z,'--',...
    out.PEDM.U2_fun(z),z,'-.',...
    UxDSV(inds),-1./(2*wavenumbers(inds)),'s',...
    U*cosd(psi),z,'k',...
    -U*sind(psi),z,'-.k','LineWidth',2);

%ylim([min(z),0]);
%ylim([-1./(2*kVals(inds(1))),0]);
ylim([-1./(2*wavenumbers(1)),0]);
xlim(1.0*[-1,1.0]);
xlabel('[m/s]');ylabel('Depth [m]');

p(3).MarkerFaceColor = 0.6*[1,1,1];
p(3).MarkerEdgeColor = 'k';
p(3).LineWidth = 1.0;

legend('PEDM_x','PEDM_y','EDM_x','U_x(z)','U_y(z)','interpreter','none','location','se');
set(gca,'fontsize',12);
drawnow;


UxE = interp1(kT,Uk(1,:)*cosd(psi),wavenumbers);
UyE = interp1(kT,-Uk(1,:)*sind(psi),wavenumbers);

outExact = find_current_depth_profile(wavenumbers(inds),UxE(inds),UyE(inds),Inf);

if saveRes
    resName = fullfile(dataDir,fileName);

    
    dsName = [locName,'/UxPEDM_3'];
    dataVar = out.PEDM.U1_fun(z);
    h5create(resName,dsName,size(dataVar));
    h5write(resName,dsName,dataVar);
        
    dsName = [locName,'/UyPEDM_3'];
    dataVar = out.PEDM.U2_fun(z);
    h5create(resName,dsName,size(dataVar));
    h5write(resName,dsName,dataVar);
    
    dsName = [locName,'/UxPEDM_exact_3'];
    dataVar = outExact.PEDM.U1_fun(z);
    h5create(resName,dsName,size(dataVar));
    h5write(resName,dsName,dataVar);
        
    dsName = [locName,'/UyPEDM_exact_3'];
    dataVar = outExact.PEDM.U2_fun(z);
    h5create(resName,dsName,size(dataVar));
    h5write(resName,dsName,dataVar);
    
end

