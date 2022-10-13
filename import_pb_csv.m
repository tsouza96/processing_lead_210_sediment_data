%% Import data from text file

cd (strcat(datafolder))
source_dir = uigetdir([]); 
myfiles = dir(source_dir);
filenames = {myfiles(:).name}';
csvfile = filenames(endsWith(filenames,'.csv'));

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 8);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["depth", "pbtotal", "upb", "ratotal", "ura", 'dbd','udbd', 'year'];
opts.VariableTypes = ['double', 'double', 'double', 'double', "double" ,'double' ,'double', 'double'];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data - CHANGE FILE NAME HERE
tbl = readtable(csvfile{1});

%% Convert to output type
depth = tbl.depth;
pbtotal = tbl.pbtotal;
upb = tbl.upb;
ratotal = tbl.ratotal;
ura = tbl.ura;
dbd = tbl.dbd;
udbd = tbl.udbd;
year = tbl.year;

%% Clear temporary variables
clear opts tbl

save('model')