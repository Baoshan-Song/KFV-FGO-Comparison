classdef state < handle
    properties
        gid      % global ID in all states        
        lid      % local ID in sliding window

        value   % state value
        status  % status: forward, backward, smooth, margin, fixed, invalid
    end
    
    methods
        function obj = state(gid, lid, value)
            obj.gid = gid;
            obj.lid = lid;            
            obj.value = value;
        end
    end
end
