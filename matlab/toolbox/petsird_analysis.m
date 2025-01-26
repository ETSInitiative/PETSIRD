function petsird_analysis(input)

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
fprintf("Number of modules of first type: %d\n", ...
    length(header.scanner.scanner_geometry.replicated_modules(1).transforms));
fprintf("Number of types of detecting elements in modules of first type: %d\n", ...
    length(header.scanner.scanner_geometry.replicated_modules(1).object.detecting_elements));
fprintf("Number of elements of first type in modules of first type: %d\n", ...
    length(header.scanner.scanner_geometry.replicated_modules(1).object.detecting_elements(1).transforms));
fprintf("Total number of 'crystals': %d\n", ...
    petsird.helpers.get_num_detecting_elements(header.scanner.scanner_geometry));
fprintf("Number of TOF bins: %d\n", header.scanner.number_of_tof_bins());
fprintf("Number of energy bins: %d\n", header.scanner.number_of_event_energy_bins());

event_energy_bin_edges = header.scanner.event_energy_bin_edges;
fprintf("Energy bin edges: " + join(repelem("%f", length(event_energy_bin_edges))) + "\n", event_energy_bin_edges);

energy_mid_points = (event_energy_bin_edges(1:end-1) + event_energy_bin_edges(2:end)) / 2;
fprintf("Event energy mid points: " + join(repelem("%f", length(energy_mid_points))) + "\n", energy_mid_points);

% sgidlut = header.scanner.detection_efficiencies.module_pair_sgidlut;
% fprintf("SGID LUT:\n");
% disp(sgidlut);

energy_1 = 0.0;
energy_2 = 0.0;
num_prompts = 0;
num_delayeds = 0;
last_time = 0;
while reader.has_time_blocks()
    time_block = reader.read_time_blocks();

    if time_block.isEventTimeBlock()
        last_time = time_block.value.start;
        num_prompts = num_prompts + length(time_block.value.prompt_events);
        if time_block.value.delayed_events ~= yardl.None
            num_delayeds = num_delayeds + length(time_block.value.delayed_events.value);
        end
        fprintf("=====================  Events in time block from %f  ==============\n", last_time);
        for event_idx = 1:length(time_block.value.prompt_events)
            event = time_block.value.prompt_events(event_idx);
            energy_1 = energy_1 + energy_mid_points(event.energy_indices(1) + 1);
            energy_2 = energy_2 + energy_mid_points(event.energy_indices(2) + 1);

            fprintf("CoincidenceEvent(detectorIds=[%d, %d], tofIdx=%d, energyIndices=[%d, %d])\n", ...
                event.detector_ids(1), event.detector_ids(2), event.tof_idx, ...
                event.energy_indices(1), event.energy_indices(2));
            mod_and_el = petsird.helpers.get_module_and_element(header.scanner.scanner_geometry, event.detector_ids);
            fprintf("    [ModuleAndElement(module=%d, el=%d), ModuleAndElement(module=%d, el=%d)]\n", ...
                mod_and_el(1).module, mod_and_el(1).el, mod_and_el(2).module, mod_and_el(2).el);
            fprintf("    efficiency: %f\n", petsird.helpers.get_detection_efficiency(header.scanner, event));
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
