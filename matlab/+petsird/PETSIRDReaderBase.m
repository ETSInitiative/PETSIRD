% This file was generated by the "yardl" tool. DO NOT EDIT.

% Definition of the stream of data
classdef PETSIRDReaderBase < handle
  properties (Access=protected)
    state_
  end

  methods
    function self = PETSIRDReaderBase()
      self.state_ = 0;
    end

    function close(self)
      self.close_();
      if self.state_ ~= 2
        expected_method = self.state_to_method_name_(self.state_);
        throw(yardl.ProtocolError("Protocol reader closed before all data was consumed. Expected call to '%s'.", expected_method));
      end
    end

    % Ordinal 0
    function value = read_header(self)
      if self.state_ ~= 0
        self.raise_unexpected_state_(0);
      end

      value = self.read_header_();
      self.state_ = 1;
    end

    % Ordinal 1
    function more = has_time_blocks(self)
      if self.state_ ~= 1
        self.raise_unexpected_state_(1);
      end

      more = self.has_time_blocks_();
      if ~more
        self.state_ = 2;
      end
    end

    function value = read_time_blocks(self)
      if self.state_ ~= 1
        self.raise_unexpected_state_(1);
      end

      value = self.read_time_blocks_();
    end

    function copy_to(self, writer)
      writer.write_header(self.read_header());
      while self.has_time_blocks()
        item = self.read_time_blocks();
        writer.write_time_blocks({item});
      end
      writer.end_time_blocks();
    end
  end

  methods (Static)
    function res = schema()
      res = petsird.PETSIRDWriterBase.schema;
    end
  end

  methods (Abstract, Access=protected)
    read_header_(self)
    has_time_blocks_(self)
    read_time_blocks_(self)

    close_(self)
  end

  methods (Access=private)
    function raise_unexpected_state_(self, actual)
      actual_method = self.state_to_method_name_(actual);
      expected_method = self.state_to_method_name_(self.state_);
      throw(yardl.ProtocolError("Expected call to '%s' but received call to '%s'.", expected_method, actual_method));
    end

    function name = state_to_method_name_(self, state)
      if state == 0
        name = "read_header";
      elseif state == 1
        name = "read_time_blocks";
      else
        name = "<unknown>";
      end
    end
  end
end
