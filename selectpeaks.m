%
% File:   selectpeaks.m
% Date:   2022-10-13
% Author: J. Himmelstein <himmelstein@unc.edu>
%
% Matlab function to extract counts for Po-209 and Po-210 based on known keV
% ranges of interest/data in CHN files. Potential for drift of peaks, thus
% re-calibration should be manually performed by plotting individual CHN
% files and re-delineating peaks
% 
% Usage: counts = selectpeaks(plotflag)
%
% Where:  
%
% plotflag  - optional argument: 1 plot Pb-210[dpm/g] vs Depth[cm]
%
% User prompted inputs will be required, follow formatting found on example
% files in repository [e.g. 2021_Detector-Background.txt with MCAnum
% Po-209	Po-209_err	Po-210	Po-210_err	Po-209_start	Po-209_end
% Po-210_start	Po-210_end as headers and 2022-09-27_WC1_Core-Data.txt with
% Depth (cm) and Fine mass added to Tube (g) as column headers]
%
% counts is a structure with length = # of samples and the following fields:
%
% depth             -  vector of depth
% counts_209_210    - raw counts within ROI(defined in det_bkgrnd_file)
% cpm_209_210       - Net counts per minute [209 210]
% Po_210_dpm_g      - Po-210 dpm/g 
% name              - Name of file / depth interval
% detector          - Detector used for sampling and ROI definition
% all               - Structure with all data produced by readchn function

% Example usage:
%
% data = selectpeaks() % will not auto-plot data (can be plotted later)
%
% or
% 
% data = selectpeaks(1) % will plot Pb-210 Profile at end of running


function counts = selectpeaks(plotflag)

if(nargin<1)		% if plotflag was not specified, then default to 0
  plotflag = 0;
end

fds = fileDatastore('*.chn', 'ReadFcn', @importdata);
fullFileNames = fds.Files;
numFiles = length(fullFileNames);

disp("Add Detector Background File (with regions of interest for each peak included)")
det_bkgrnd_file = uigetfile('/*.TXT;*.TSV'); %add detector background file
det_bkgrnds = readtable(det_bkgrnd_file);
disp("Add Core Data File [with Depth(cm), Sediment added to vial (g)]")
core_data_file = uigetfile('/*.TXT;*.TSV');
core_data = readtable(core_data_file);

% Loop over all files in d reading them in

for i = 1 : numFiles

%     fprintf('Now reading file %s\n', fullFileNames{i}); %Uncomment if you
%     wish to see files being read in (for Quality Control)
    
    mcadat = readchn(fullFileNames{i},1);

%     plot(mcadat.energy,mcadat.count);   % plot counts versus energy
%     axis('tight'); ylabel('Decay Counts'); xlabel('energy [KeV]'); grid on;
%     title(mcadat.filename);
%     zoom on;
%     movegui('west');

    det_idx= find(ismember(det_bkgrnds.MCAnum,(mcadat.mcanum),'rows'));
    idx_209 = (mcadat.keV>det_bkgrnds.Po_209_start(det_idx) & mcadat.keV<det_bkgrnds.Po_209_end(det_idx));
    idx_210 = (mcadat.keV>det_bkgrnds.Po_210_start(det_idx) & mcadat.keV<det_bkgrnds.Po_210_end(det_idx));

    counts_209_210 = [sum(mcadat.count(idx_209)),sum(mcadat.count(idx_210))];

    data.depth = core_data.Depth_cm_(i);
    data.counts_209_210 = counts_209_210;
    data.cpm_209_210 = (counts_209_210/mcadat.livetime)*60 - [det_bkgrnds.Po_209(det_idx) det_bkgrnds.Po_210(det_idx)];
    data.Po_210_dpm_g = ((data.cpm_209_210(2)/data.cpm_209_210(1))*(10*1))/(core_data.FineMassAddedToTube_g_(i)); % (Po210/Po209)*(tracer_activity_dpm*pipette_mL)/(samplemass_mg/1000)
    data.cpm_error_209_210 = [det_bkgrnds.Po_209_err(det_idx) det_bkgrnds.Po_210_err(det_idx)];
    %     sqrt(((Po_209_net_cpm/Po_210_net_cpm)*sqrt((([Po_210_net_cpm_error/Po_210_net_cpm])^2)+(([Po_09_net_cpm_error/Po_209_net_cpm])^2))^2/[Po-210_dpm/g_]^2)+(([tracer_activity_dpm_error^2)/([tracer_activity_dpm]^2))+(([pipette(mL)error]^2)/([pipette(mL)]^2]))+(((.005*[sample_mass_mg_error])^2)/([sample_mass_mg(i)]^2)))
    data.Po_210_dpm_g_error = sqrt(((data.cpm_209_210(1)/data.cpm_209_210(2))*sqrt(((data.cpm_error_209_210(2)/data.cpm_209_210(2))^2)+((data.cpm_error_209_210(1)/data.cpm_209_210(1))^2))^2/data.Po_210_dpm_g(1)^2)+((.719859^2)/(10^2))+((.001812^2)/(1^2))+(((.005*core_data.FineMassAddedToTube_g_(i))^2)/((core_data.FineMassAddedToTube_g_(i))^2)));
    %     data.yield = =((Po_209_net_cpm)/(det_bkgrnds.total_alpha_3_8(det_idx)))/(tracer_activity*pipette_vol))*100);
    data.yield = (((data.cpm_209_210(1)/(det_bkgrnds.total_alpha_3_8(det_idx)))/(10*1))*100);
    data.name = extractAfter((string(mcadat.filename)),106); % change character # to extract after depending on your filename length
    data.detector = mcadat.detector;
    data.all = mcadat;
    
    % for given iteration add data to count structure
    counts(i) = data;
end

if(plotflag == 1)
% If Desired, Plot Po-210 (dpm/g) data (y-values need changing)
    clf;
    plot([counts.Po_210_dpm_g],[counts.depth],'o');
    set(gca, 'YDir','reverse');
    title('Pb-210 Profile');
    axis('tight'); ylabel('Depth [cm]'); xlabel('Pb-210 [dpm/g]'); grid on;
end
