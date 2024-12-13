% This file was generated by the "yardl" tool. DO NOT EDIT.

classdef RigidTransformationSerializer < yardl.binary.RecordSerializer
  methods
    function self = RigidTransformationSerializer()
      field_serializers{1} = yardl.binary.FixedNDArraySerializer(yardl.binary.Float32Serializer, [4, 3]);
      self@yardl.binary.RecordSerializer('petsird.RigidTransformation', field_serializers);
    end

    function write(self, outstream, value)
      arguments
        self
        outstream (1,1) yardl.binary.CodedOutputStream
        value (1,1) petsird.RigidTransformation
      end
      self.write_(outstream, value.matrix);
    end

    function value = read(self, instream)
      fields = self.read_(instream);
      value = petsird.RigidTransformation(matrix=fields{1});
    end
  end
end
