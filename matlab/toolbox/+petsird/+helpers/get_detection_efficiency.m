function eff = get_detection_efficiency(scanner, event)
% Compute the detection efficiency for a coincidence event

arguments
    scanner (1,1) petsird.ScannerInformation
    event (1,1) petsird.CoincidenceEvent
end

eff = 1;

% Per det_el efficiencies
det_el_efficiencies = scanner.detection_efficiencies.det_el_efficiencies;
if det_el_efficiencies ~= yardl.None
    eff = eff .* (det_el_efficiencies(event.energy_indices(1)+1, event.detector_ids(1)+1) .* ...
            det_el_efficiencies(event.energy_indices(2)+1, event.detector_ids(2)+1));
end

% Per module-pair efficiencies
module_pair_efficiencies_vector = scanner.detection_efficiencies.module_pair_efficiencies_vector;
if module_pair_efficiencies_vector ~= yardl.None
    module_pair_SGID_LUT = scanner.detection_efficiencies.module_pair_sgidlut;
    assert(module_pair_SGID_LUT ~= yardl.None);
    mod_and_els = petsird.helpers.get_module_and_element(scanner.scanner_geometry, event.detector_ids);
    assert(length(scanner.scanner_geometry.replicated_modules) == 1);
    SGID = module_pair_SGID_LUT(mod_and_els(1).module+1, mod_and_els(2).module+1);
    assert(SGID >= 0);
    module_pair_efficiencies = module_pair_efficiencies_vector(SGID+1);
    assert(module_pair_efficiencies.sgid == SGID);
    eff = eff * module_pair_efficiencies.values( ...
        event.energy_indices(2)+1, ...
        mod_and_els(2).el+1, ...
        event.energy_indices(1)+1, ...
        mod_and_els(1).el+1 ...
    );
end

end
