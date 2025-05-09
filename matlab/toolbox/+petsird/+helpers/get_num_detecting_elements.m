function num_det_els = get_num_detecting_elements(scanner_geometry)
% Compute total number of detecting elements in the scanner

arguments
    scanner_geometry (1,1) petsird.ScannerGeometry
end

num_det_els = 0;
for rep_mod_idx = 1:length(scanner_geometry.replicated_modules)
    rep_module = scanner_geometry.replicated_modules(rep_mod_idx);
    det_els = rep_module.object.detecting_elements;
    num_det_els = num_det_els + length(det_els.transforms) .* length(rep_module.transforms);
end

end
