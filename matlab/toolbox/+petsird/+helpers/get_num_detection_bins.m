function num_bins = get_num_detection_bins(scanner, type_of_module)
% Compute total number of detecting bins in a spcific module-type in the scanner

arguments
    scanner (1,1) petsird.ScannerInformation
    type_of_module (1,1) petsird.TypeOfModule
end

rep_module = scanner.scanner_geometry.replicated_modules(type_of_module);
det_els = rep_module.object.detecting_elements;
energy_bin_edges = scanner.event_energy_bin_edges(type_of_module);
num_bins = length(det_els.transforms) .* length(rep_module.transforms) ...
  * energy_bin_edges.number_of_bins();

end
