function num_det_els = get_num_detecting_elements(scanner_geometry, type_of_module)
% Compute total number of detecting elements in a spcific module-type in the scanner

arguments
    scanner_geometry (1,1) petsird.ScannerGeometry
    type_of_module (1,1) petsird.TypeOfModule
end

rep_module = scanner_geometry.replicated_modules(type_of_module);
det_els = rep_module.object.detecting_elements;
num_det_els = length(det_els.transforms) .* length(rep_module.transforms);

end
