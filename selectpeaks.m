%
% File:   selectpeaks.m
% Date:   2023-02-22
% Author: J. Himmelstein <himmelstein@unc.edu>
%
% Matlab function to extract counts for Po-209 and Po-210 based on known keV
% ranges of interest/data in CHN files. Potential for drift of peaks, thus
% re-calibration should be manually performed by plotting individual CHN
% files and re-delineating peaks
%
% Usage: counts = selectpeaks(filetype, plotflag)
%
% Where:
%
% filetype - required input: "chn" or "spe" - based on your filetype (.chn is
% maestro native filetype, .spe is an ASCII file format)
% plotflag  - optional argument: 1 = plot Pb-210[dpm/g] vs Depth[cm]
%
% User prompted inputs will be required, follow formatting found on example
% files in repository [e.g. 2021_Detector-Background.txt with MCAnum
% Po-209	Po-209_err	Po-210	Po-210_err as headers and 2022-09-27_WC1_Core-Data.txt with
% Depth (cm) and Fine mass added to Tube (g) as column headers]
%
% counts is a structure with length = # of samples and the following fields:
%
% depth             - vector of depth
% counts_209_210    - raw counts within ROI(defined in det_bkgrnd_file)
% cpm_209_210       - Net counts per minute [209 210]
% Po_210_dpm_g      - Po-210 dpm/g
% name              - Name of file / depth interval
% detector          - Detector used for sampling and ROI definition
% all               - Structure with all data produced by readchn function

% Example usage:
%
% data = selectpeaks("spe") % if your files from the detector are saved as .spe
% and because no plotflag speficifed, will not auto-plot data (can be plotted later)
%
% or
%
% data = selectpeaks("chn",1) % if your files from the detector are saved
% as .chn and will plot Pb-210 Profile at end of running


function counts = selectpeaks(filetype, plotflag)

if(nargin<2)		% if plotflag was not specified, then default to 0
  plotflag = 0;
end

if any((filetype ~= "chn") & (filetype ~= "spe"))
    fprintf("ERROR: Filetype not supported or too few inputs")
    return
end

if filetype == "chn"
    disp("Select folder with .CHN files to count peaks from (one core only)")
    spectrum_files = uigetdir(); %add CHN spectrum file directory
    file_list = dir(fullfile(spectrum_files,'/*.chn'));

elseif filetype == "spe"
    disp("Select folder with .SPE files to count peaks from (one core only)")
    spectrum_files = uigetdir(); %add SPE spectrum file directory
    file_list = dir(fullfile(spectrum_files,'/*.spe'));

end

file_table = struct2table(file_list);
fullFileNames = append(file_table.folder,'\',file_table.name);
fullFileNames = string(fullFileNames);
numFiles = height(fullFileNames);

disp("Add Detector Background File")
det_bkgrnd_file = uigetfile('/*.TXT;*.CSV;*.TSV'); %add detector background file
% det_bkgrnd_file = 'C:\Users\Rodriguez Lab\Desktop\Josh-Desktop\Marsh_Levee_Project\Wards_Creek\2011_Detector-Backgrounds.txt';
det_bkgrnds = readtable(det_bkgrnd_file);

disp("Add Core Data File [with Depth(cm), Sediment added to vial (g)]")
core_data_file = uigetfile('/*.TXT;*.TSV');

% folder = 'C:\Users\Rodriguez Lab\Desktop\Josh-Desktop\Marsh_Levee_Project\Wards_Creek';
% core_data_file = fullfile(folder, '2022-09-27_WC-1_Core-Data.csv');

core_data = readtable(core_data_file);

% Loop over all files in d reading them in

for i = 1 : numFiles

%     fprintf('Now reading file %s\n', fullFileNames{i}); %Uncomment if you
%     wish to see files being read in (for Quality Control)
    if filetype == "chn"
        mcadat = readchn(fullFileNames(i),1);
    elseif filetype =="spe"
        mcadat = readspe(fullFileNames(i),1);
    end

%     plot(mcadat.energy,mcadat.count);   % plot counts versus energy
%     axis('tight'); ylabel('Decay Counts'); xlabel('energy [KeV]'); grid on;
%     title(mcadat.filename);
%     zoom on;
%     movegui('west');


    det_idx= find(ismember(det_bkgrnds.MCAnum,(mcadat.mcanum),'rows'));

    % Find the index of the first element greater than 5 (Start of 209
    % Peak)
    index1 = find(mcadat.count > 15, 1, 'first');

    % Find the index of the first element less than 1 after the first
    % element greater than 5 (End of 209 Peak)
    index2 = find(mcadat.count(index1+1:end) < 1, 1, 'first') + index1;

    % Find the index of the first element greater than 5 after index 2
    % (Start of 210 peak)
    index3 = find(mcadat.count(index2:end) > 15, 1, 'first') + index2;

    % Find the index of the first element less than 1 after the third index
    % found above (End of 210 Peak)
    index4 = find(mcadat.count(index3:end) < 1, 1, 'first') + index3;

    counts_209_210 = [sum(mcadat.count(index1:index2)),sum(mcadat.count(index3:index4))];

    % Assign variables to data structure
    data.depth = core_data.Depth_cm_(i); % depth data (must be at least as many core depths as samples processed, else error)
    data.counts_209_210 = counts_209_210; % 209 and 210 counts
    data.detector = mcadat.detector; % detector data

    % Calculate CPM and CPM error for 209/210
    cpm_209_210 = (counts_209_210/mcadat.livetime)*60 - [det_bkgrnds.Po_209(det_idx) det_bkgrnds.Po_210(det_idx)];
    cpm_error_209_210 = [sqrt((cpm_209_210(1)/mcadat.livetime/60)+(det_bkgrnds.Po_209_err(det_idx)/mcadat.livetime/60)+(det_bkgrnds.Po_209_err(det_idx)/mcadat.livetime/60)) sqrt((cpm_209_210(2)/mcadat.livetime/60)+(det_bkgrnds.Po_210_err(det_idx)/mcadat.livetime/60)+(det_bkgrnds.Po_210_err(det_idx)/mcadat.livetime/60))];
    data.cpm_209_210 = cpm_209_210; % assign cpm to data structure
    data.cpm_error_209_210 = cpm_error_209_210; % assign cpm error to data structure

    % Calculate Po-210 DPM/g
    Po_210_dpm_g = ((cpm_209_210(2)/cpm_209_210(1))*(10*1))/(core_data.FineMassAddedToTube_g_(i));
    data.Po_210_dpm_g = Po_210_dpm_g; % assign Po-210 DPM/g to data structure

    % Calculate Po-210 DPM/g error
    Po_210_dpm_g_error = sqrt(((cpm_209_210(1)/cpm_209_210(2))*sqrt(((cpm_error_209_210(2)/cpm_209_210(2))^2)+((cpm_error_209_210(1)/cpm_209_210(1))^2))^2/Po_210_dpm_g(1)^2)+((.719859^2)/(10^2))+((.001812^2)/(1^2))+(((.005*core_data.FineMassAddedToTube_g_(i))^2)/((core_data.FineMassAddedToTube_g_(i))^2)));
    data.Po_210_dpm_g_error = Po_210_dpm_g_error; % assign Po-210 DPM/g error to data structure

    % Calculate yield
    yield = (((cpm_209_210(1)/(det_bkgrnds.total_alpha_3_8(det_idx)))/(10*1))*100);
    data.yield = yield; % assign yield to data structure

    % Assign filename to data structure
    [filepath,name,ext] = fileparts(mcadat.filename);
    data.name = name;
    
    % Assign all data from mcadat to data structure
    data.all = mcadat;
   
    % for given iteration add data to count structure
    counts(i) = data;
end

if(plotflag == 1)
% If Desired, Plot Po-210 (dpm/g) data (y-values need changing)
    clf;
    plot([counts.Po_210_dpm_g],[counts.depth],'o--');
    set(gca, 'YDir','reverse');
    title('Pb-210 Profile');
    axis('tight'); ylabel('Depth [cm]'); xlabel('Pb-210 [dpm/g]'); grid on;
end

% % Uncomment to write your data to a table that can be used in age-depth modeling, otherwise work
% % with matlab structure created from this function
outputdata = table([counts.depth]', [counts.Po_210_dpm_g]', [counts.Po_210_dpm_g_error]');
outputdata.Properties.VariableNames = ["depth", "pbtotal", "upb"];
writetable(outputdata,'counts.csv')


