
% PROP_PMA  Generates the prop struct used to summarize CPMA parameters.
%  
%  PROP = prop_pma() creates a default PMA properties structure for use in 
%  evaluating transfer function. This is equiavlent to prop_pma('olfert').
%  
%  PROP = prop_pma(SPEC) add a string specifying parameter set. 
%  
%  AUTHOR: Timothy Sipkens, 2019-06-26

function [prop] = prop_pma(spec)

% if type of property specification is not given, use 'olfert'.
if ~exist('spec', 'var'); spec = []; end
if isempty(spec); spec = 'olfert'; end

%-- Default mass-mobility information -------------%
% Default, Dm = 3 corresponds to spheres (no change in effective density).
prop.Dm = 3; % mass-mobility exponent
prop.m0 = 4.7124e-25; % mass-mobility pre-factor

% Universal relation (Olfert and Rogak)
% prop.Dm = 2.48;
% prop.m0 = 2.9280e-24;


% Read case-specific properties from YAML file.
% Uses read_yaml(...) method given as a subfunction in this file.
[fd, ~] = fileparts(mfilename('fullpath'));  % get current folder
prop = read_yaml([fd, filesep, 'prop', filesep, spec, '.yaml'], prop);


%-- Parameters related to CPMA geometry ----------------------------------%
prop.rc = (prop.r2 + prop.r1) / 2;
prop.r_hat = prop.r1 / prop.r2;
prop.del = (prop.r2 - prop.r1) / 2; % half gap width

prop.A = pi * (prop.r2 ^ 2 - prop.r1 ^ 2); % cross sectional area of APM
prop.v_bar = prop.Q / prop.A; % average flow velocity


%-- For diffusion --------------------------------------------------------%
kB = 1.3806488e-23; % Boltzmann's constant
prop.D = @(B) kB .* prop.T .* B; % diffusion coefficient


% Fill mass-mobility relation equivalents.
prop = prop_massmob(prop);

end




% READ_YAML  Simple utility to read simple YAML files.
%  
%  AUTHOR: Lloyd Russell, 2017
%  REPO: https://github.com/llerussell/ReadYAML
%  MODIFIED: Timothy Sipkens, 2021-04-19
function prop = read_yaml(file_path, prop)

% Parse inputs (currently never invoked).
if ~exist('prop', 'var'); prop = []; end
if isempty(prop); prop = struct(); end

% If file of presets does not exist.
if ~isfile(file_path)
    error(['PMA properties not available for provided case. ',...
            'Try a different string in the `prop_pma(str)` call.']);
end


% Read file line by line.
fid = fopen(file_path, 'r');
data = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');
fclose(fid);

% Remove empty lines.
data = deblank(data{1});
data(cellfun('isempty', data)) = [];

% Prepare final results structure.
results = [];

% Parse the contents (line by line).
for i = 1:numel(data)
    
    % extract this line
    thisLine = data{i};
    
    % ignore if this line is a comment
    if strcmpi(thisLine(1), '#')
        continue
    end
    
    % find the seperator between key and value
    sepIndex = find(thisLine==':', 1, 'first');
    
    % get the key name (remove whitespace)
    key = strtrim(thisLine(1:sepIndex-1));  
    
    % get the value, ignoring any comments (remove whitespace)
    value = strsplit(thisLine(sepIndex+1:end), '#');
    value = strtrim(value{1}); 
    
    % attempt to convert value to numeric type
    [convertedValue, success] = str2num(value);
    if success
        value = convertedValue;
    end
    
    % store the key and value in the results
    results.(key) = value;
end

% Append to supplied structure.
% Overwrite duplicate fields. 
f_results = fields(results);
for ii = 1:length(f_results)
    prop.(f_results{ii}) = results.(f_results{ii});
end

end


