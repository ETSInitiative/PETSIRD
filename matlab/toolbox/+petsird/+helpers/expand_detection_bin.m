function expanded_detection_bin = expand_detection_bin(scanner_geometry, type_of_module, detection_bin)
% Find ExpandedDetectionBin for a detection_bin

arguments
    scanner_geometry (1,1) petsird.ScannerGeometry
    type_of_module (1,1) petsird.TypeOfModule
    detection_bin (1,1) petsird.DetectionBin
end

rep_module = scanner_geometry.replicated_modules(type_of_module);
num_el_per_module = length(rep_module.object.detecting_elements.transforms);
energy_bin_edges = rep_module.event_energy_bin_edges(type_of_module);
num_en = energy_bin_edges.number_of_bins();
expanded_detection_bin = petsird.ExpandedDetectionBin('moduleIndex', idivide(detection_bin, num_el_per_module * num_en), ...
                                                      'elementIndex', mod(idivide(detection_bin, num_en), num_el_per_module), ...
                                                      'energyIndex', mod(detection_bin, num_en));

end
