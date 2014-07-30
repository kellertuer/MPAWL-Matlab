function setDebugLevel(varargin)
% setDebugLevel(type,level)
%   Set the global debug level of type to the value level, where
%   the type can be ommited, which yields setting types 'all' and 'any'
%      'all' is the upper bound, 'any' the lower bound for all types.
%   
%   INPUT
%       type  : (optional) a string of the type of debug the leve is set
%       level : level the type is set to
%
% ---
% LMDLO ~ R. Bergmann ~ 2014-04-21, last edit: 2014-04-22

    global LMDLO_struct;
    % if debug was not active before, calling setDebugLevel activates it.
    initializeDebug();
    if (nargin<1) || (nargin>2)
        warning('Wrong number of setDebugLevel (',num2str(nargin),...
            '). Ingoring function call.');
        return;
    end
    if (nargin == 1)
        setDebugLevel('all',varargin{1});
        setDebugLevel('any',varargin{1});
        return;
    end
    type = varargin{1};
    level = varargin{2};
    if ~any(strcmp(type,LMDLO_struct.types))
        %type not yet included
        newind = length(LMDLO_struct.types)+1;
        LMDLO_struct.types{newind} = type;
        LMDLO_struct.levels(newind) = level;
        if strcmp(type,'any') && level > getDebugLevel('all')
            warning(['Trying to set ''any'' (',num2str(level),') > ''all'' (',num2str(getDebugLevel('all')),') is not possible, using the value of ''all''']);
            LMDLO_struct.levels(newind) = getDebugLevel('all');
        end
        debug('text',3,'Text',['Debug Level of ',type,' set to ',num2str(LMDLO_struct.levels(newind))]);
    else
        %update type
        thisind = find(strcmp(type,LMDLO_struct.types),1);
        LMDLO_struct.levels(thisind) = level;
        LMDLO_struct.levels(thisind) = level;
        if strcmp(type,'any') && level > getDebugLevel('all')
            warning(['Trying to set ''any'' (',num2str(level),') > ''all'' (',num2str(getDebugLevel('all')),') is not possible, using the value of ''all''']);
            LMDLO_struct.levels(thisind) = getDebugLevel('all');
        end
        debug('text',3,'Text',['Debug Level of ',type,' set to ',num2str(LMDLO_struct.levels(thisind))]);
    end
end