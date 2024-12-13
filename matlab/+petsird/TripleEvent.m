% This file was generated by the "yardl" tool. DO NOT EDIT.

classdef TripleEvent < handle
  % All information about a triple event specified as identifiers or indices (i.e. discretized).
  properties
    detector_ids
    % timing differences (converted to mm) w.r.t. first event, stored as
    % indices into the tofBinEdges field in the ScannerInformation
    % Note: only 2, corresponding to the arrival time differences of the second and third detectorId
    % listed w.r.t. the first detectorId
    tof_indices
    % indices for each single event into the energyBinEdges field in the ScannerInformation
    energy_indices
  end

  methods
    function self = TripleEvent(kwargs)
      arguments
        kwargs.detector_ids = repelem(uint32(0), 3);
        kwargs.tof_indices = repelem(uint32(0), 2);
        kwargs.energy_indices = repelem(uint32(0), 3);
      end
      self.detector_ids = kwargs.detector_ids;
      self.tof_indices = kwargs.tof_indices;
      self.energy_indices = kwargs.energy_indices;
    end

    function res = eq(self, other)
      res = ...
        isa(other, "petsird.TripleEvent") && ...
        isequal({self.detector_ids}, {other.detector_ids}) && ...
        isequal({self.tof_indices}, {other.tof_indices}) && ...
        isequal({self.energy_indices}, {other.energy_indices});
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
      elem = petsird.TripleEvent();
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
