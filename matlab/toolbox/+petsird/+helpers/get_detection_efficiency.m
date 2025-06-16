function eff = get_detection_efficiency(scanner, type_of_module_pair, event)
% Compute the detection efficiency for a coincidence event

arguments
    scanner (1,1) petsird.ScannerInformation
    type_of_module_pair (1,1) petsird.TypeOfModule
    event (1,1) petsird.CoincidenceEvent
end

eff = 1;

type_of_module_1 = type_of_module_pair(1)
type_of_module_2 = type_of_module_pair(2)

% Per detection_bin efficiencies
detection_bin_efficiencies = scanner.detection_efficiencies.detection_bin_efficiencies;
if detection_bin_efficiencies ~= yardl.None
    detection_bin_efficiencies_1 = ...
        detection_efficiencies.detection_bin_efficiencies(type_of_module_pair_1)
    detection_bin_efficiencies_2 = ...
        detection_efficiencies.detection_bin_efficiencies(type_of_module_pair_2)

  eff = eff .* (detection_bin_efficiencies_1(event.detection_bins(1)+1) .* ...
            detection_bin_efficiencies_2(event.detection_bins(2)+1));
end

% Per module-pair efficiencies
module_pair_efficiencies_vector = scanner.detection_efficiencies.module_pair_efficiencies_vector;
if module_pair_efficiencies_vector ~= yardl.None
    module_pair_SGID_LUT = scanner.detection_efficiencies.module_pair_sgidlut;
    assert(module_pair_SGID_LUT ~= yardl.None);
    expanded_det_bin_1 = petsird.helpers.expand_detection_bin(scanner.scanner_geometry, type_of_module_1, event.detection_bins(1));
    expanded_det_bin_2 = petsird.helpers.expand_detection_bin(scanner.scanner_geometry, type_of_module_2, event.detection_bins(2));
    SGID = module_pair_SGID_LUT(type_of_module_1)(type_of_module_2)(expanded_det_bin_1.module+1, expanded_det_bin_2.module+1);
    if SGID < 0
        eff = 0;
        return
    end
    module_pair_efficiencies = module_pair_efficiencies_vector(type_of_module_1)(type_of_module_2)(SGID+1);
    assert(module_pair_efficiencies.sgid == SGID);
    eff = eff * module_pair_efficiencies.values( ...
        event.detection_bins(2).energy_idx+1, ...
        expanded_det_bin_2.el+1, ...
        event.detection_bins(1).energy_idx+1, ...
        expanded_det_bin_1.el+1 ...
    );
end

end
