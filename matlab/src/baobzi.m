
classdef baobzi < handle
    properties (SetAccess = private, Hidden = true)
        objectHandle; % Handle to the underlying C++ class instance
    end
    methods
        function this = baobzi(varargin)
            this.objectHandle = baobzi_mex('new', varargin{:});
        end

        function delete(this)
            baobzi_mex('free', this.objectHandle);
        end

        function varargout = eval(this, varargin)
            [varargout{1:nargout}] = baobzi_mex('eval', this.objectHandle, varargin{:});
        end

        function varargout = save(this, varargin)
            [varargout{1:nargout}] = baobzi_mex('save', this.objectHandle, varargin{:});
        end

        function varargout = restore(this, varargin)
          [varargout{1:nargout}] = baobzi_mex('restore', this.objectHandle, varargin{:});
        end
    end
end
