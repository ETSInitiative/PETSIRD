function expanded_det_bin = expand_detection_bin(scanner_geometry, type_of_module, detection_bin)
% Find ExpandedDetectionBin for a detection_bin

arguments
    scanner_geometry (1,1) petsird.ScannerGeometry
    type_of_module (1,1) petsird.TypeOfModule
    detection_bin (1,1) petsird.DetectionBin
end

rep_module = scanner_geometry.replicated_modules(type_of_module);
num_el_per_module = length(rep_module1.object.detecting_elements.transforms);

expanded_det_bin = struct('module', idivide(detection_bin.det_el_idx, num_el_per_module), 'el', mod(detection_bin.det_el_idx, num_el_per_module));

end
