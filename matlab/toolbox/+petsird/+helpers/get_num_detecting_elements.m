function num_det_els = get_num_detecting_elements(scanner, type_of_module)
% Compute total number of detecting elements in a specific module-type in the scanner

arguments
    scanner (1,1) petsird.Scanner
    type_of_module (1,1) petsird.TypeOfModule
end

rep_module = scanner.scanner_geometry.replicated_modules(type_of_module);
det_els = rep_module.object.detecting_elements;
num_det_els = length(det_els.transforms) .* length(rep_module.transforms);

end
