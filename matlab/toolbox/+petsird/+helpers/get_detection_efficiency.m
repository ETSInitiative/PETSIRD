function eff = get_detection_efficiency(scanner, event)
% Compute the detection efficiency for a coincidence event

arguments
    scanner (1,1) petsird.ScannerInformation
    event (1,1) petsird.CoincidenceEvent
end

eff = 1;

% Per detection_bin efficiencies
detection_bin_efficiencies = scanner.detection_efficiencies.detection_bin_efficiencies;
if detection_bin_efficiencies ~= yardl.None
    eff = eff .* (detection_bin_efficiencies(event.detection_bins(1).energy_idx+1, event.detection_bins(1).det_el_idx+1) .* ...
            detection_bin_efficiencies(event.detection_bins(2).energy_idx+1, event.detection_bins(2).det_el_idx+1));
end

% Per module-pair efficiencies
module_pair_efficiencies_vector = scanner.detection_efficiencies.module_pair_efficiencies_vector;
if module_pair_efficiencies_vector ~= yardl.None
    module_pair_SGID_LUT = scanner.detection_efficiencies.module_pair_sgidlut;
    assert(module_pair_SGID_LUT ~= yardl.None);
    mod_and_els = petsird.helpers.get_module_and_element(scanner.scanner_geometry, event.detection_bins);
    assert(length(scanner.scanner_geometry.replicated_modules) == 1);
    SGID = module_pair_SGID_LUT(mod_and_els(1).module+1, mod_and_els(2).module+1);
    if SGID < 0
        eff = 0;
        return
    end
    module_pair_efficiencies = module_pair_efficiencies_vector(SGID+1);
    assert(module_pair_efficiencies.sgid == SGID);
    eff = eff * module_pair_efficiencies.values( ...
        event.detection_bins(2).energy_idx+1, ...
        mod_and_els(2).el+1, ...
        event.detection_bins(1).energy_idx+1, ...
        mod_and_els(1).el+1 ...
    );
end

end
