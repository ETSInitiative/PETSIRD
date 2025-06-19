function detection_bin = make_detection_bin(scanner_geometry, type_of_module, detection_bin)
% Find DetectionBin for a ExpandedDetectionBin

arguments
    scanner_geometry (1,1) petsird.ScannerGeometry
    type_of_module (1,1) petsird.TypeOfModule
    expanded_detection_bin (1,1) petsird.ExpandedDetectionBin
end

rep_module = scanner_geometry.replicated_modules(type_of_module);
num_el_per_module = length(rep_module.object.detecting_elements.transforms);
energy_bin_edges = rep_module.event_energy_bin_edges(type_of_module);
num_en = energy_bin_edges.number_of_bins();
% use shorter name
bin = expanded_detection_bin;
detection_bin = petsird.DetectionBin(bin.energy_index + (bin.element_index + bin.module_index * num_el_per_module) * num_en);

end
