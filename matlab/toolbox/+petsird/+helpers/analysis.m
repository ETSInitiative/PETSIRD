function analysis(input)

arguments
    input (1,1) string
end

reader = petsird.binary.PETSIRDReader(input);

header = reader.read_header();

if header.exam ~= yardl.None
    fprintf("Subject ID: %s\n", header.exam.subject.id);
end

fprintf("Scanner name: %s\n", header.scanner.model_name);

fprintf("Types of modules: %d\n", length(header.scanner.scanner_geometry.replicated_modules));
type_of_module = 1;
fprintf("Number of modules of first type: %d\n", ...
    length(header.scanner.scanner_geometry.replicated_modules(type_of_module).transforms));
fprintf("Number of elements in modules of first type: %d\n", ...
    length(header.scanner.scanner_geometry.replicated_modules(type_of_module).object.detecting_elements.transforms));
fprintf("Total number of 'crystals': %d\n", ...
    petsird.helpers.get_num_detecting_elements(header.scanner.scanner_geometry));
tof_bin_edges = header.scanner.tof_bin_edges[type_of_module][type_of_module];
num_tof_bins = tof_bin_edges.number_of_bins()

event_energy_bin_edges = header.scanner.event_energy_bin_edges(type_of_module);
num_event_energy_bins = event_energy_bin_edges.number_of_bins();
event_energy_bin_edges = event_energy_bin_edges.edges;
fprintf("Number of TOF bins: %d\n", num_tof_bins);
fprintf("Number of energy bins: %d\n", num_event_energy_bins);

fprintf("Energy bin edges: " + join(repelem("%f", length(event_energy_bin_edges))) + "\n", event_energy_bin_edges);

energy_mid_points = (event_energy_bin_edges(1:end-1) + event_energy_bin_edges(2:end)) / 2;
fprintf("Event energy mid points: " + join(repelem("%f", length(energy_mid_points))) + "\n", energy_mid_points);

sgidlut = header.scanner.detection_efficiencies.module_pair_sgidlut;
fprintf("SGID LUT:\n");
disp(sgidlut);

energy_1 = 0.0;
energy_2 = 0.0;
num_prompts = 0;
num_delayeds = 0;
last_time = 0;
while reader.has_time_blocks()
    time_block = reader.read_time_blocks();

    if time_block.isEventTimeBlock()
        % TODO just doing one module-type ATM
        type_of_module_pair = petsird.TypeOfModulePair((1, 1))
        type_of_module_1 = type_of_module_pair(1);
        type_of_module_2 = type_of_module_pair(2);
        last_time = time_block.value.start;
        num_prompts = num_prompts + length(time_block.value.prompt_events(type_of_module_1)(type_of_module_2);
        if time_block.value.delayed_events ~= yardl.None
            num_delayeds = num_delayeds + length(time_block.value.delayed_events(type_of_module_1)(type_of_module_2).value);
        end
        fprintf("=====================  Events in time block from %f  ==============\n", last_time);
        for event_idx = 1:length(time_block.value.prompt_events(type_of_module_1)(type_of_module_2))
            event = time_block.value.prompt_events(type_of_module_1)(type_of_module_2)(event_idx);
            energy_1 = energy_1 + energy_mid_points(event.detection_bins(1).energy_idx + 1);
            energy_2 = energy_2 + energy_mid_points(event.detection_bins(2).energy_idx + 1);

            fprintf("CoincidenceEvent(detElIndices=[%d, %d], tofIdx=%d, energyIndices=[%d, %d])\n", ...
                event.detection_bins(1).det_el_idx, event.detection_bins(2).det_el_idx, event.tof_idx, ...
                event.detection_bins(1).energy_idx,  event.detection_bins(2).energy_idx);
            expanded_det_bin_1 = petsird.helpers.expand_detection_bins(header.scanner.scanner_geometry, type_of_module_1, event.detection_bins(1))
            expanded_det_bin_2 = petsird.helpers.expand_detection_bins(header.scanner.scanner_geometry, type_of_module_2, event.detection_bins(2))
            fprintf("    [ExpandedDetectionBin(module=%d, el=%d), ExpandedDetectionBin(module=%d, el=%d)]\n", ...
                expanded_det_bin_1.module, expanded_det_bin_1.el, expanded_det_bin_2.module, expanded_det_bin_2.el);
            fprintf("    efficiency: %f\n", petsird.helpers.get_detection_efficiency(header.scanner, type_of_module_pair, event));
        end
    end
end

fprintf("Last time block at %f ms\n", last_time);
fprintf("Number of prompt events: %d\n", num_prompts);
fprintf("Number of delayed events: %d\n", num_delayeds);
if num_prompts > 0
    fprintf("Average energy_1: %f\n", energy_1 / num_prompts)
    fprintf("Average energy_2: %f\n", energy_2 / num_prompts)
end

reader.close();

end
