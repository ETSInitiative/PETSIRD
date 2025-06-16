function expanded_det_bins = expand_detection_bins(scanner_geometry, type_of_module, detection_bins)
% Find ExpandedDetectionBin for a list of detection_bins

arguments
    scanner_geometry (1,1) petsird.ScannerGeometry
    type_of_module (1,1) petsird.TypeOfModule
    detection_bins (:,1) petsird.DetectionBin
end

rep_module = scanner_geometry.replicated_modules(type_of_module);
num_el_per_module = length(rep_module1.object.detecting_elements.transforms);

expanded_det_bins = arrayfun(@(id) struct('module', idivide(id.det_el_idx, num_el_per_module), 'el', mod(id.det_el_idx, num_el_per_module)), detection_bins);

end
