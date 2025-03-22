function mod_and_el = get_module_and_element(scanner_geometry, detection_bins)
% Find ModuleAndElement for a list of detection_bins

arguments
    scanner_geometry (1,1) petsird.ScannerGeometry
    detection_bins (:,1) petsird.DetectionBin
end

assert(length(scanner_geometry.replicated_modules) == 1);
rep_module = scanner_geometry.replicated_modules(1);
assert(length(rep_module.object.detecting_elements) == 1);
num_el_per_module = length(rep_module.object.detecting_elements(1).ids);

mod_and_el = arrayfun(@(id) struct('module', idivide(id.det_el_idx, num_el_per_module), 'el', mod(id.det_el_idx, num_el_per_module)), detection_bins);

end
