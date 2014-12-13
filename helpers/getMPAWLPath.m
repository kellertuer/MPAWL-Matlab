function folder = getMPAWLPath()
% getMVDPPath()
%    returns the base path of the MVDP Toolbox and returns an error
%    message, if the Toolbox is not initialized yet.
%
% ---
% R. Bergmann ~ 2014-12-13

folder = fileparts(which('initMPAWL.m'));
assert(~isempty(folder),...
    'MPAWL not found in path, please run initMPAWL.m first.');
end

