function [state] = result2state(results,time_ind)
%UNTITLED Summary of this function goes here
%   state.u
%   state.ux
%   state.uy
%   state.uz
%   state.time

    switch class(results)
        case {'pde.StationaryResults','hhxt.HHXTResults','hhxt.HHXTStationaryResults'}
        state.time = 0;
        state.u = results.NodalSolution';
        state.ux = results.XGradients';
        state.uy = results.YGradients';
        if ~results.IsTwoD
            state.uz = results.ZGradients';
        end
        case {'pde.TransientResults','hhxt.HHXTTransientResults'}
        state.time = results.SolutionTimes(time_ind);
        [cols,rows] = size(results.NodalSolution(:,:,1));
        state.u = zeros(rows,cols,length(time_ind));
        state.ux = state.u; state.uy = state.u;
        if ~results.IsTwoD; state.uz = state.u; end
        
        for t = 1:length(time_ind)
            state.u(:,:,t) = results.NodalSolution(:,:,t)';
            state.ux(:,:,t) = results.XGradients(:,:,t)';
            state.uy(:,:,t) = results.YGradients(:,:,t)';
            if ~results.IsTwoD
                state.uz(:,:,t) = results.ZGradients(:,:,t)';
            end
        end
        
        otherwise
            error('unrecognized result')
    end



end

