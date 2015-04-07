%NESTEDFOR actPos for a for-loop with dynamic depth
%   I = NESTEDFOR(start,end) constructs the actPos. The index vectors produced by the actPos will fulfill
%   start <= index <= end componentwise
%
%   I = I.next() advances the index by one step. The last component of the index is running the fastest.
%
%   [I,index] = I.next() advances the index by one step and returns the resulting index.
%
%   index = I.actPos gives the current index.
%
%   I.finished() returns 0 if index <= end componentwise and 1 otherwise
%
% ---
% MPAWL; D. Merkert ~ 2014-08-14, last edit: 2014-08-30 (R. Bergmann)

%   EXAMPLE:
%   start = [1 2];
%   end   = [3 5];
%   I = nestedFor(start,end);
%   while(I.hasNext())
%     display(I.next());
%   end

classdef nestedFor < handle %inherit from handle a nicer way to work with obj
    properties
        startV
        endV
        actPos
    end
    methods
        function obj=nestedFor(startV_,endV_)
            assert(isrow(startV_),'start should be a row vector');
            assert(isrow(endV_),'end should be a row vector');
            assert(size(startV_,2) == size(endV_,2), 'start and end should have the same dimensions');
            assert(all(startV_ <= endV_), 'It should hold start(i) <= end(i) for all i');
            obj.startV = startV_;
            obj.endV = endV_;
            obj.actPos = NaN;
        end
        function it = next(obj)
            if isnan(obj.actPos)
                obj.actPos = obj.startV;
            else
                if ~obj.hasNext()
                    error('No next iterate existing.');
                end
                %search for smallest entry less then end to increment
                nextstep = find(obj.actPos<obj.endV,1);
                obj.actPos(nextstep) = obj.actPos(nextstep)+1;
                % reset all previous ones that reached end
                obj.actPos(1:(nextstep-1)) = obj.startV(1:(nextstep-1)); 
                if any( (obj.actPos>obj.endV)|(obj.actPos<obj.startV))
                    error('Error during iteration, actPos invalid. Did you miss to check hasNext?');
                end
            end
            it = obj.actPos;
        end
        function res = hasNext(obj)
            res = any(obj.actPos < obj.endV) | isnan(obj.actPos);
        end
    end
end