function lvl = getDebugLevel( type )
%   getDebugLevel(type)
%       get Debug Level of a certain type, 0 if not set
%   INPUT
%       type    : type a of debug
%
%   OUTPUT
%       lvl     : level in {1,2,3} of that debug type, if set,
%                   this level is bounded by `all` (if set) from above
%                   this lebel is bounded by `and` (if set) from below
%                 0 is returned, if none of these three is set
% ---
% LMDLO, R. Bergmann, 2014-04-21
    global LMDLO_active LMDLO_struct;
    if (isempty(LMDLO_active)) %not initailized by any setDebugLevel
        warning('getDebugLevel() should bot be called without initialization first');
        lvl = -1;
        return;
    end
    lvl=0;
    if any(strcmp(type,LMDLO_struct.types))
        %For this type exists an explicit value
        lvl = LMDLO_struct.levels(find(strcmp(type,LMDLO_struct.types),1));
        all_lvl=lvl;
        any_lvl=lvl;
    else
        all_lvl=Inf; any_lvl=0;
    end
    if any(strcmp('all',LMDLO_struct.types)) %is all set? remember that level
        all_lvl = LMDLO_struct.levels(find(strcmp('all',LMDLO_struct.types),1));
    end
    if any(strcmp('any',LMDLO_struct.types)) %is any set? remember that level
        any_lvl = LMDLO_struct.levels(find(strcmp('any',LMDLO_struct.types),1));
    end
    lvl = min(max(any_lvl,lvl),all_lvl);
end