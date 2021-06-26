function vndatastatic = importfile(filename, dataLines)
%IMPORTFILE Import data from a text file
%  VNDATASTATIC = IMPORTFILE(FILENAME) reads data from text file
%  FILENAME for the default selection.  Returns the data as a table.
%
%  VNDATASTATIC = IMPORTFILE(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  vndatastatic = importfile("E:\Studium\M2-Inertialnavigation\Uebung\2\vn-data-static.csv", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 2021-06-11 22:35:51

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 59);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["GpsCycle", "GpsWeek", "GpsToW", "TimeStartup", "TimeSyncIn", "SyncInCnt", "UnCompMagX", "UnCompMagY", "UnCompMagZ", "UnCompAccX", "UnCompAccY", "UnCompAccZ", "UnCompGyroX", "UnCompGyroY", "UnCompGyroZ", "Temperature", "Pressure", "DeltaTime", "DeltaThetaX", "DeltaThetaY", "DeltaThetaZ", "DeltaVelX", "DeltaVelY", "DeltaVelZ", "MagX", "MagY", "MagZ", "AccX", "AccY", "AccZ", "GyroX", "GyroY", "GyroZ", "AhrsStatus", "Yaw", "Pitch", "Roll", "Quat0", "Quat1", "Quat2", "Quat3", "MagN", "MagE", "MagD", "AccN", "AccE", "AccD", "LinAccX", "LinAccY", "LinAccZ", "LinAccN", "LinAccE", "LinAccD", "YawU", "PitchU", "RollU", "YawRate", "PitchRate", "RollRate"];
opts.VariableTypes = ["string", "string", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["GpsCycle", "GpsWeek", "GpsToW", "YawRate", "PitchRate", "RollRate"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["GpsCycle", "GpsWeek", "GpsToW", "YawRate", "PitchRate", "RollRate"], "EmptyFieldRule", "auto");

% Import the data
vndatastatic = readtable(filename, opts);

end