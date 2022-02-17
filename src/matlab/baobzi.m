
classdef baobzi < handle
  properties (SetAccess = private, Hidden = true)
    objectHandle; % Handle to the underlying C++ class instance
  end
  methods
    function this = baobzi(varargin)
      nargs = length(varargin);
      if nargs < 1
        error("Baobzi needs at least one argument: 'restore' or 'new'.")
      end

      if strcmp(varargin{1}, 'new')
        this.objectHandle = baobzi_mex('new', varargin{2:end});
      elseif strcmp(varargin{1}, 'restore')
        this.objectHandle = baobzi_mex('restore', varargin{2:end});
      else
        error("Baobzi's first argument should be either 'restore', or 'new'.")
      end
    end

    function delete(this)
      baobzi_mex('free', this.objectHandle);
    end

    function varargout = eval(this, varargin)
      [varargout{1:nargout}] = baobzi_mex('eval', this.objectHandle, varargin{:});
    end

    function varargout = stats(this, varargin)
      [varargout{1:nargout}] = baobzi_mex('stats', this.objectHandle, varargin{:});
    end

    function varargout = save(this, varargin)
      [varargout{1:nargout}] = baobzi_mex('save', this.objectHandle, varargin{:});
    end
  end
end
