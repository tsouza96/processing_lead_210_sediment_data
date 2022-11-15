%
% File:   readspe.m
% Date:   2022-11-15
% Author: J. Himmelstein <himmelstein@unc.edu>, reference to  I. Chuang <ichuang@mit.edu>
%
% Matlab function to read Alpha counter data saved in the general
%  "spe" format.
%
% Usage:  mcadat = readspe(filename,plotflag)
%
% Where:  
%
% filename  - name of ASCII format data file (foo.spe) to read in
% plotflag  - optional argument: 1 plot count vs chan, 2 plot count vs energy
% mcadat    - output data structure
%
% mcadat is a structure with the following fields:
%
% chan      - vector of channel numbers
% count     - vector of counts in each channel number
% energy    - vector of energies (given by the software energy calib) [KeV]
% mcanum    - MCA number (1-4)
% realtime  - acquisition real-time duration [sec]
% livetime  - acquisition live-time duration [sec]
% dtstamp   - date-time stamp (date,month,year,hour,minutes,sec) as a string
% startchan - channel offset to starting channel
% nchan     - number of channels
% econv     - energy conversion factor [KeV]/[chan]
% detector  - string describing the detector
% filename  - filename which data was stored in
% 
% Example usage:
%
% mcadat = readspe('2022-09-15_Det1_WC-1_0-1cm.spe')
% figure(1)
% plot(mcadat.chan,mcadat.count)     % plot counts versus channel number
% axis('tight'); ylabel('counts'); xlabel('channel'); grid on
% title(mcadat.filename)
% figure(2)
% plot(mcadat.energy,mcadat.count)   % plot counts versus energy
% axis('tight'); ylabel('counts'); xlabel('energy [KeV]'); grid on
%
% or
% 
% readspe('2022-09-15_Det1_WC-1_0-1cm.spe', 1)

function mcadat = readspe(fname, plotflag)

if(nargin<2)		% if plotflag was not specified, then default to 0
  plotflag = 0;
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 7);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["SPEC_ID", "VarName2", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7"];
opts.VariableTypes = ["string", "string", "string", "string", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% Specify variable properties
opts = setvaropts(opts, ["SPEC_ID", "VarName2", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["SPEC_ID", "VarName2", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7"], "EmptyFieldRule", "auto");

% Import the data
filedata= readmatrix(fname, opts);

%% Clear temporary variables
clear opts

%% Begin parsing file into structure

mcadat.filename = fname;
mcadat.mcanum = str2double(filedata(4,2));
mcadat.realtime = str2num(filedata(10,1));
mcadat.livetime = str2num(filedata(10,2));
mcadat.dtstamp = filedata(8,1) +" " + filedata(8,2);
mcadat.startchan = str2double(filedata(12,1));
mcadat.nchan = str2double(filedata(12,2))+1;
mcadat.count = str2double(filedata(12:mcadat.nchan+11,1));	% counts
%filemark = fread(fp,1,'int16');		% file marker
mcadat.econv = str2num(filedata(end-3,2));		% energy conversion factor
mcadat.detector = filedata(5,2)+" "+filedata(5,3)+" "+filedata(5,4);
mcadat.detector_long =filedata(5,2)+" "+filedata(5,3)+" "+filedata(5,4)+" "+filedata(5,5)+" "+filedata(5,6)+" "+filedata(5,7);
% make vector of channels
mcadat.chan = mcadat.startchan:mcadat.startchan+mcadat.nchan-1; 
% make vector of keV [JDH on 2022-10-07]
mcadat.keV = [3.10:.01:44.05];
% make vector of energies
mcadat.energy = mcadat.chan * mcadat.econv;

% plot, if requested

if(plotflag==1)
  clf;
  plot(mcadat.chan,mcadat.count);	% plot counts versus channel
  axis('tight'); ylabel('counts'); xlabel('channel'); grid on;
  title(mcadat.filename);
  zoom on;
end

if(plotflag==2)
  clf;
  plot(mcadat.energy,mcadat.count);   % plot counts versus energy
  axis('tight'); ylabel('counts'); xlabel('energy [KeV]'); grid on;
  title(mcadat.filename);
  zoom on;
end
  
