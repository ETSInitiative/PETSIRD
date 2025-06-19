function detection_bins = make_detection_bins(scanner, type_of_module, expanded_detection_bins)
% Find DetectionBin for a list of ExpandedDetectionBins

arguments
    scanner_geometry (1,1) petsird.ScannerInformation
    type_of_module (1,1) petsird.TypeOfModule
    expanded_detection_bins (:,1) petsird.ExpandedDetectionBin
end

% TODO slow but avoids re-implementation
detection_bins = arrayfun(@(bin) make_detection_bin(scanner, type_of_module, bin)), expanded_detection_bins);

end
