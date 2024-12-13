% This file was generated by the "yardl" tool. DO NOT EDIT.

classdef TimeFrameInformation < handle
  % A sequence of time intervals (could be consecutive)
  properties
    time_frames
  end

  methods
    function self = TimeFrameInformation(kwargs)
      arguments
        kwargs.time_frames = petsird.TimeInterval.empty();
      end
      self.time_frames = kwargs.time_frames;
    end

    function res = number_of_time_frames(self)
      res = length(self.time_frames);
      return
    end


    function res = eq(self, other)
      res = ...
        isa(other, "petsird.TimeFrameInformation") && ...
        isequal({self.time_frames}, {other.time_frames});
    end

    function res = ne(self, other)
      res = ~self.eq(other);
    end

    function res = isequal(self, other)
      res = all(eq(self, other));
    end
  end

  methods (Static)
    function z = zeros(varargin)
      elem = petsird.TimeFrameInformation();
      if nargin == 0
        z = elem;
        return;
      end
      sz = [varargin{:}];
      if isscalar(sz)
        sz = [sz, sz];
      end
      z = reshape(repelem(elem, prod(sz)), sz);
    end
  end
end
