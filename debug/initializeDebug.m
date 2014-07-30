function wasReset = initializeDebug(reset)
% initializeDebug(reset) check whether Debug is initialized, and initialize
% if necessary. The optional parameter reset can be used to reset the Debug
% to its initial state, which is the same upon first activation.
%
%  INPUT:
%       reset : (optional) if set, it resets the LMDLO to its initial
%       working state
%
%  OUTPUT:
%        wasReset : if the Debug System was (re)initialized, the function
%        returns true, else false
%
% LMDLO ~ R. Bergmann ~ 2014-04-21

    global LMDLO_active LMDLO_struct;
    if (... %if one of the global struct elements is missing
        (~isfield(LMDLO_struct,'types'))...
        || (~isfield(LMDLO_struct,'levels'))... %or we should reset
        || ((nargin > 0) && (reset))...
    )
    % Initialize debug
        LMDLO_active = true;
        LMDLO_struct.types = {};
        LMDLO_struct.levels(1) = 0;
        wasReset = true;   
    else
        wasReset = false;
end