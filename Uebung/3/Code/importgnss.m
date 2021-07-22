function vn310gnss = importgnss(filename, dataLines)
%IMPORTFILE Import data from a text file
%  VN310GNSS = IMPORTFILE(FILENAME) reads data from text file FILENAME
%  for the default selection.  Returns the data as a table.
%
%  VN310GNSS = IMPORTFILE(FILE, DATALINES) reads data for the specified
%  row interval(s) of text file FILENAME. Specify DATALINES as a
%  positive scalar integer or a N-by-2 array of positive scalar integers
%  for dis-contiguous row intervals.
%
%  Example:
%  vn310gnss = importfile("E:\Studium\M2-Inertialnavigation\Uebung\3\Code\vn310-gnss.csv", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 2021-07-22 13:54:37

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 137);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["GpsCycle", "GpsWeek", "GpsTow", "TimeTimeStartup", "TimeGpsTow", "TimeGpsWeek", "TimeTimeUtcyear", "TimeTimeUtcmonth", "TimeTimeUtcday", "TimeTimeUtchour", "TimeTimeUtcmin", "TimeTimeUtcsec", "TimeTimeUtcms", "TimeTimeStatustimeOk", "TimeTimeStatusdateOk", "TimeTimeStatusutcTimeValid", "GNSS1UTCyear", "GNSS1UTCmonth", "GNSS1UTCday", "GNSS1UTChour", "GNSS1UTCmin", "GNSS1UTCsec", "GNSS1UTCms", "GNSS1Tow", "GNSS1Week", "GNSS1NumSats", "GNSS1Fix", "GNSS1PosLlalatitude", "GNSS1PosLlalongitude", "GNSS1PosLlaaltitude", "GNSS1PosEcefX", "GNSS1PosEcefY", "GNSS1PosEcefZ", "GNSS1VelNedN", "GNSS1VelNedE", "GNSS1VelNedD", "GNSS1VelEcefX", "GNSS1VelEcefY", "GNSS1VelEcefZ", "GNSS1PosUN", "GNSS1PosUE", "GNSS1PosUD", "GNSS1VelU", "GNSS1TimeU", "GNSS1TimeInfoStatustimeOk", "GNSS1TimeInfoStatusdateOk", "GNSS1TimeInfoStatusutcTimeValid", "GNSS1TimeInfoLeapSeconds", "GNSS1gDOP", "GNSS1pDOP", "GNSS1tDOP", "GNSS1vDOP", "GNSS1hDOP", "GNSS1nDOP", "GNSS1eDOP", "AttYawPitchRollY", "AttYawPitchRollP", "AttYawPitchRollR", "AttQuaternionw", "AttQuaternionx", "AttQuaterniony", "AttQuaternionz", "AttYprUY", "AttYprUP", "AttYprUR", "INSInsStatusMode", "INSInsStatusGpsFix", "INSInsStatusErrorIMU", "INSInsStatusErrorMagPres", "INSInsStatusErrorGNSS", "INSInsStatusGpsHeadingIns", "INSInsStatusGpsCompass", "INSPosLlalatitude", "INSPosLlalongitude", "INSPosLlaaltitude", "INSPosEcefX", "INSPosEcefY", "INSPosEcefZ", "INSVelBodyX", "INSVelBodyY", "INSVelBodyZ", "INSVelNedN", "INSVelNedE", "INSVelNedD", "INSVelEcefX", "INSVelEcefY", "INSVelEcefZ", "INSMagEcefX", "INSMagEcefY", "INSMagEcefZ", "INSAccelEcefX", "INSAccelEcefY", "INSAccelEcefZ", "INSLinearAccelEcef", "INSLinearAccelEcefY", "INSLinearAccelEcefZ", "INSPosU", "INSVelU", "GNSS2UTCyear", "GNSS2UTCmonth", "GNSS2UTCday", "GNSS2UTChour", "GNSS2UTCmin", "GNSS2UTCsec", "GNSS2UTCms", "GNSS2Tow", "GNSS2Week", "GNSS2NumSats", "GNSS2Fix", "GNSS2PosLlalatitude", "GNSS2PosLlalongitude", "GNSS2PosLlaaltitude", "GNSS2PosEcefX", "GNSS2PosEcefY", "GNSS2PosEcefZ", "GNSS2VelNedN", "GNSS2VelNedE", "GNSS2VelNedD", "GNSS2VelEcefX", "GNSS2VelEcefY", "GNSS2VelEcefZ", "GNSS2PosUN", "GNSS2PosUE", "GNSS2PosUD", "GNSS2VelU", "GNSS2TimeU", "GNSS2TimeInfoStatustimeOk", "GNSS2TimeInfoStatusdateOk", "GNSS2TimeInfoStatusutcTimeValid", "GNSS2TimeInfoLeapSeconds", "GNSS2gDOP", "GNSS2pDOP", "GNSS2tDOP", "GNSS2vDOP", "GNSS2hDOP", "GNSS2nDOP", "GNSS2eDOP"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
vn310gnss = readtable(filename, opts);

end