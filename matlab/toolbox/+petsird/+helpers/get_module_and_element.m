function mod_and_el = get_module_and_element(scanner_geometry, scanner_det_ids)
% Find ModuleAndElement for a list of detector_ids

arguments
    scanner_geometry (1,1) petsird.ScannerGeometry
    scanner_det_ids (:,1) uint32
end

assert(length(scanner_geometry.replicated_modules) == 1);
rep_module = scanner_geometry.replicated_modules(1);
assert(length(rep_module.object.detecting_elements) == 1);
num_el_per_module = length(rep_module.object.detecting_elements(1).ids);

mod_and_el = arrayfun(@(id) struct('module', idivide(id, num_el_per_module), 'el', mod(id, num_el_per_module)), scanner_det_ids);

end
