function expanded_detection_bins = expand_detection_bins(scanner, type_of_module, detection_bins)
% Find ExpandedDetectionBin for a list of detection_bins

arguments
    scanner_geometry (1,1) petsird.ScannerInformation
    type_of_module (1,1) petsird.TypeOfModule
    detection_bins (:,1) petsird.DetectionBin
end

% TODO slow but avoids re-implementation
expanded_detection_bins = arrayfun(@(bin) expand_detection_bin(scanner, type_of_module, bin)), detection_bins);

end
