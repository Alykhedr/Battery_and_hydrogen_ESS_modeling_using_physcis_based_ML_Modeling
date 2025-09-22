function [time, P_el, P_fc] = yearly_power_profile(csvFile, varargin)
% YEARLY_LOAD_PROFILE  Read PV/FC profiles from a tidy CSV and return model-ready vectors.
%   [time, P_el, P_fc] = yearly_load_profile('profiles_PV1MW_FCscaled_60min.csv', ...)
%
% Input CSV columns (kW):  time(ISO-8601 UTC), P_pv, P_load
%
% Options (name/value):
%   'Start'     : datetime or char, inclusive (default: first timestamp in file)
%   'End'       : datetime or char, inclusive (default: last timestamp in file)
%   'TimeZone'  : 'UTC' | 'Europe/Berlin' | etc. (default: 'UTC')
%   'AsHours'   : true → return hours since Start; false → return datetime vector (default: true)
%
% Outputs:
%   time  : column vector; hours since Start (double) or datetime (tz-aware)
%   P_el  : PV electric power [kW]   (nonnegative, column)
%   P_fc  : FC demand power [kW]     (nonnegative, column)
%
% Example:
%   [t, Pel, Pfc] = yearly_load_profile('profiles_PV1MW_FCscaled_60min.csv', ...
%                       'Start','2019-01-01T00:00:00Z', 'End','2019-12-31T23:00:00Z', ...
%                       'TimeZone','Europe/Berlin', 'AsHours', true);

  % ---- defaults ----
  if nargin < 1 || isempty(csvFile)
    csvFile = 'profiles_PV1MW_FCscaled_60min.csv';
  end
  p.Start    = [];
  p.End      = [];
  p.TimeZone = 'UTC';
  p.AsHours  = true;

  % parse name/value args
  if ~isempty(varargin)
    for k = 1:2:numel(varargin)
      name = varargin{k};
      val  = varargin{k+1};
      switch lower(name)
        case 'start',    p.Start = val;
        case 'end',      p.End   = val;
        case 'timezone', p.TimeZone = val;
        case 'ashours',  p.AsHours  = logical(val);
        otherwise, error('Unknown option: %s', name);
      end
    end
  end

  % ---- read table ----
  opts = detectImportOptions(csvFile, 'NumHeaderLines', 0);
  must = {'time','P_pv','P_load'};
  assert(all(ismember(must, opts.VariableNames)), 'CSV missing required columns time/P_pv/P_load.');
  T = readtable(csvFile, opts);

  % ---- parse time (stored as UTC strings) ----
  tsUTC = string(T.time);
  tUTC  = datetime(tsUTC, 'InputFormat',"yyyy-MM-dd'T'HH:mm:ss'Z'", 'TimeZone','UTC');

  % optional window
  if ~isempty(p.Start)
    if ~isduration(p.Start) && ~isdatetime(p.Start)
      p.Start = datetime(char(p.Start), 'InputFormat',"yyyy-MM-dd'T'HH:mm:ss'Z'", 'TimeZone','UTC');
    elseif isdatetime(p.Start) && isempty(p.Start.TimeZone)
      p.Start.TimeZone = 'UTC';
    end
  else
    p.Start = tUTC(1);
  end

  if ~isempty(p.End)
    if ~isduration(p.End) && ~isdatetime(p.End)
      p.End = datetime(char(p.End), 'InputFormat',"yyyy-MM-dd'T'HH:mm:ss'Z'", 'TimeZone','UTC');
    elseif isdatetime(p.End) && isempty(p.End.TimeZone)
      p.End.TimeZone = 'UTC';
    end
  else
    p.End = tUTC(end);
  end

  mask = (tUTC >= p.Start) & (tUTC <= p.End);
  assert(any(mask), 'Selected window has no rows.');
  tUTC  = tUTC(mask);
  Ppv   = T.P_pv(mask);
  Pload = T.P_load(mask);

  % hygiene
  Ppv   = max(0, Ppv(:));
  Pload = max(0, Pload(:));

  % ---- timezone convert for output clock ----
  tOut = tUTC;
  if ~strcmpi(p.TimeZone, 'UTC')
    tOut.TimeZone = p.TimeZone;
  end

  % ---- outputs aligned with your model naming ----
  P_el = Ppv;    % PV electric power to electrolyzer side [kW]
  P_fc = Pload;  % FC demand [kW]

  % ---- time as hours-since-start or datetime ----
  if p.AsHours
    time = hours(tUTC - tUTC(1));
    time = time(:);
  else
    time = tOut(:);
  end

  % sanity checks
  assert(isequal(numel(time), numel(P_el), numel(P_fc)), 'Size mismatch in outputs.');
end
