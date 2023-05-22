% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Omid Sani, Maryam Shanechi
%   Shanechi Lab, University of Southern California, 2017
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dFilter Data buffer
%   Usage example:
%       flt = dFilter;

classdef dFilter < handle & matlab.mixin.Copyable
    
    properties (SetAccess = private)
        b        % numerator
        a        % denumerator
        state    % state of the filter
        causal   % if true will aplly causally
        Fs       % Intended sampling rate of data for the filter
    end
    
    properties 
        
    end
    
    events
        
    end
    
    methods
        function obj = dFilter(b, a, causal, Fs)
            % Initializes a filter
            % Inputs: 
            % - (1) b: numerator [or a structure of the inputs, 'a', 'b', ...]
            % - (2) a: denumerator 
            % - (3) causal (default: true). If true, will apply filter
            %           causally and store state
            % - (4) Fs (default: []). Intended data Fs. Just for future
            %           reference
            % Outputs: 
            % - (none)
            if nargin < 1, b = 1; end
            if nargin < 2, a = 1; end
            if nargin < 3, causal = true; end
            if nargin < 4, Fs = []; end
            if isstruct(b)
                settings = b;
                if ~isfield(settings, 'b'), settings.b = 1; end
                if isfield(settings, 'b'), b = settings.b; end
                if isfield(settings, 'a'), a = settings.a; end
                if isfield(settings, 'causal'), causal = settings.causal; end
                if isfield(settings, 'Fs'), Fs = settings.Fs; end
            end
            obj.b = b;
            obj.a = a;
            obj.setCausal( causal );
            obj.state = [];
            obj.setFs( Fs );
        end
        function [data] = apply( obj, data )
            % Applies filter to data
            % Inputs: 
            % (1) x: filter will be applied to each column separately
            % Outputs: 
            % (1) data: filtered data
            
            if ~isempty(obj.state)&&(obj.getStateDim()~=size(data, 2))
                obj.reset(); % Reset if new data dim doesn't match state
            end
            isNaN = isnan(data);
            isNaNMin = min(isNaN, [], 2); 
            isNaNMax = max(isNaN, [], 2); 
            hasSameNaNsInDims = isequal(isNaNMin, isNaNMax);
            if obj.causal
                if hasSameNaNsInDims 
                    hasNoNaNs = ~any(isNaN, 2);
                    [data(hasNoNaNs, :), obj.state] = filter(obj.b, obj.a, data(hasNoNaNs, :), obj.state);
                else
                    for ci = 1:size(data, 2)
                        if isempty(obj.state)
                            thisState = [];
                        else
                            thisState = obj.state(:, ci);
                        end
                        [data(~isNaN(:, ci), ci), thisState] = filter(obj.b, obj.a, data(~isNaN(:, ci), ci), thisState);
                        if isempty(obj.state)
                            obj.state = zeros(numel(thisState), size(data, 2));
                        end
                        obj.state(:, ci) = thisState;
                    end
                end
            else
                if hasSameNaNsInDims % All chans have NaNs in the same samples
                    hasNoNaNs = ~any(isNaN, 2);
                    data(hasNoNaNs, :) = filtfilt(obj.b, obj.a, data(hasNoNaNs, :));
                else
                    for ci = 1:size(data, 2)
                        data(~isNaN(:, ci), ci) = filtfilt(obj.b, obj.a, data(~isNaN(:, ci), ci));
                    end
                end
            end
        end
        function setCausal( obj, causal )
            % Sets the causality of the filter
            % Inputs: 
            % - (1) causal: true/false will set causality to this
            % Outputs: 
            % - none
            obj.causal = causal;
        end
        function setFs( obj, fs )
            % Sets the intended data sampling rate of the filter
            % Inputs: 
            % - (1) Fs: intended data sampling rate
            % Outputs: 
            % - none
            obj.Fs = fs;
        end
        function reset( obj )
            % Resets the filter stats
            % Inputs: 
            % - none
            % Outputs: 
            % - none
            obj.state = [];
        end
        function d = getStateDim( obj )
            % Returns the dimension of data in the current filter state
            % Inputs: 
            % - none
            % Outputs: 
            % - d: dimension of data in current fiter state. 0 if no
            %       state is stored. 
            d = size(obj.state, 2);
        end
        function h = fvtool( obj )
            % Opens fvtool GUI for inspecting filter characteristics
            % Inputs: 
            % - none
            % Outputs: 
            % - (1) h: handle to the fvtool object. Can be used to change
            %           'FrequencyScale' to 'log', etc
            h = fvtool(obj.b, obj.a);
            h.Fs = obj.Fs;
        end
    end
    
    methods(Static)
        
    end
    
end

