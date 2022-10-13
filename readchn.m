
% File:   readchn.m
% Date:   29-Aug-03, 2022-10-13
% Author: I. Chuang <ichuang@mit.edu>, edits by J. Himmelstein <himmelstein@unc.edu>
%
% Matlab function to read Maestro MCA data which was saved in the
% native "chn" format.
%
% Usage:  mcadat = readchn(filename,plotflag)
%
% Where:  
%
% filename  - name of binary-format MCA data file (foo.chn) to read in
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
% mcadat = readchn('2003-08-26-compton1-45-deg-scatter-ok.chn')
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
% readchn('2003-08-26-compton1-45-deg-scatter-ok.chn',1)

function mcadat = readchn(fname,plotflag)

if(nargin<2)		% if plotflag was not specified, then default to 0
  plotflag = 0;
end

fp = fopen(fname,'r');
mcadat.filename = fname;
mcadat.filetype = fread(fp,1,'int16');
mcadat.mcanum = fread(fp,1,'int16');
mcadat.segment = fread(fp,1,'int16');
mcadat.dtstamp = '                         ';
mcadat.dtstamp(16:17) = fread(fp,2,'uchar');	% sec
mcadat.realtime = fread(fp,1,'int32')/50;
mcadat.livetime = fread(fp,1,'int32')/50;
mcadat.dtstamp(1:2) = fread(fp,2,'uchar');	% date
mcadat.dtstamp(3:5) = fread(fp,3,'uchar');	% month
mcadat.dtstamp(6:8) = fread(fp,3,'uchar');	% year
mcadat.dtstamp(10:11) = fread(fp,2,'uchar');	% hour
mcadat.dtstamp(13:14) = fread(fp,2,'uchar');	% minutes
mcadat.dtstamp(12) = ':';
mcadat.dtstamp(15) = ':';
mcadat.startchan = fread(fp,1,'int16');
mcadat.nchan = fread(fp,1,'int16');
mcadat.count = fread(fp,mcadat.nchan,'int32');	% counts
%filemark = fread(fp,1,'int16');		% file marker
dum = fread(fp,1,'float32');			% dummy data
dum = fread(fp,1,'float32');			% dummy data
mcadat.econv = fread(fp,1,'float32');		% energy conversion factor
buf = fread(fp,500,'uchar');			% read rest of buffer
mcadat.detector = char(buf(256-10:256-10+32)');

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
  
