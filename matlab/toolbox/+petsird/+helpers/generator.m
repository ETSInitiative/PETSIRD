generator(output)

arguments
    output (1,1) string
end

% These are constants for now
config.NUMBER_OF_EVENT_ENERGY_BINS = 3;
config.NUMBER_OF_TOF_BINS = 300;
config.RADIUS = 400;
config.CRYSTAL_LENGTH = [20, 4, 4];
config.NUM_CRYSTALS_PER_MODULE = [2, 4, 7];
config.NUM_MODULES_ALONG_RING = 20;
config.NUM_MODULES_ALONG_AXIS = 2;
config.MODULE_AXIS_SPACING = (config.NUM_CRYSTALS_PER_MODULE(3) + 4) * config.CRYSTAL_LENGTH(3);
config.NUMBER_OF_TIME_BLOCKS = 6;
config.COUNT_RATE = 500; % 1/ms
config.EVENT_TIME_BLOCK_DURATION = 1; % ms

writer = petsird.binary.PETSIRDWriter(output);

header = get_header(config);
writer.write_header(header);

for t = 0:config.NUMBER_OF_TIME_BLOCKS-1
    time_interval = petsird.TimeInterval(
                start=t * config.EVENT_TIME_BLOCK_DURATION,
                stop=(t + 1) * config.EVENT_TIME_BLOCK_DURATION);
    average_num = config.EVENT_TIME_BLOCK_DURATION * config.COUNT_RATE;
    % NOTE: Need Statistics and Machine Learning Toolbox for Poisson distribution functions
    % num_prompts_this_block = poissrnd(average_num);
    num_prompts_this_block = randi(config.average_num);
    type_of_module_pair = petsird.TypeOfModulePair((1, 1));
    prompts_this_block = [[get_events(header, type_of_module_pair, num_prompts_this_block, config)]];

    tb = petsird.TimeBlock.EventTimeBlock(...
        petsird.EventTimeBlock(time_interval=time_interval, prompt_events=prompts_this_block) ...
    );
    writer.write_time_blocks(tb)
end

writer.end_time_blocks();
writer.close();

end

function h = get_header(cfg)
    subject = petsird.Subject(id="123456");
    institution = petsird.Institution(name="University of Somewhere", address="Didcot, Oxfordshire, Ox11 0DE, UK");
    h = petsird.Header();
    h.exam = petsird.ExamInformation(subject=subject, institution=institution);
    h.scanner = get_scanner_info(cfg);
end

function scanner = get_scanner_info(cfg)
    scanner_geometry = get_scanner_geometry(cfg);

    num_types_of_modules = scanner_geometry.number_of_module_types()
    % TODO This code does not yet take multiple module-types into account
    assert(num_types_of_modules == 1);

    % TOF info (in mm)
    tofBinEdges = petsird.BinEdges(edges=single(linspace(-cfg.RADIUS, cfg.RADIUS, cfg.NUMBER_OF_TOF_BINS + 1)));
    energyBinEdges = petsird.BinEdges(edges=single(linspace(430, 650, cfg.NUMBER_OF_EVENT_ENERGY_BINS + 1)));
    % In this example, there is only 1 module-type, therefore we need 1x1 "nested lists"
    allTofBinEdges = [[tofBinEdges]]
    tofResolution = [[9.4]]  % in mm
    % 1-element list for energy info
    allEnergyBinEdges = [energyBinEdges]
    energyResolutionAt511 = [0.11]  % as fraction of 511

    % We need energy bin info before being able to construct the detection
    % efficiencies, so we first construct a scanner without the efficiencies
    scanner = petsird.ScannerInformation(...
        model_name="PETSIRD_TEST", ...
        scanner_geometry=scanner_geometry, ...
        tof_bin_edges=tofBinEdges, ...
        tof_resolution=tofResolution, ...
        event_energy_bin_edges=energyBinEdges, ...
        energy_resolution_at_511=energyResolutionAt511
    );

    % Now add the efficiencies
    scanner.detection_efficiencies = get_detection_efficiencies(scanner, cfg);

    scanner_info.coincidence_policy = petsird.CoincidencePolicy.REJECT_MULTIPLES;
    scanner_info.delayed_coincidences_are_stored = false;
    scanner_info.triple_events_are_stored = false;

end

function geometry = get_scanner_geometry(cfg)
    detector_module = get_detector_module(cfg);

    rep_module = petsird.ReplicatedDetectorModule(object=detector_module);
    for i = 0:cfg.NUM_MODULES_ALONG_RING-1
        angle = 2 * pi * i / cfg.NUM_MODULES_ALONG_RING;
        for ax_mod = 0:cfg.NUM_MODULES_ALONG_AXIS-1
            transform = petsird.RigidTransformation(matrix=...
                transpose([...
                    cos(angle), -sin(angle), 0, 0; ...
                    sin(angle), cos(angle), 0, 0; ...
                    0, 0, 1, cfg.MODULE_AXIS_SPACING * ax_mod; ...
                ]) ...
            );
            rep_module.transforms = [rep_module.transforms transform];
        end
    end

    geometry = petsird.ScannerGeometry(replicated_modules=[rep_module]);
end

function detector = get_detector_module(cfg)
    % Return a module of NUM_CRYSTALS_PER_MODULE cuboids
    crystal = get_crystal(cfg);
    rep_volume = petsird.ReplicatedBoxSolidVolume(object=crystal);
    N0 = cfg.NUM_CRYSTALS_PER_MODULE(1);
    N1 = cfg.NUM_CRYSTALS_PER_MODULE(2);
    N2 = cfg.NUM_CRYSTALS_PER_MODULE(3);
    for rep0 = 0:N0-1
        for rep1 = 0:N1-1
            for rep2 = 0:N2-1
                transform = petsird.RigidTransformation(matrix=...
                    transpose([...
                        1, 0, 0, cfg.RADIUS + rep0 * cfg.CRYSTAL_LENGTH(1); ...
                        0, 1, 0, (rep1 - N1 / 2) * cfg.CRYSTAL_LENGTH(2); ...
                        0, 0, 1, (rep2 - N2 / 2) * cfg.CRYSTAL_LENGTH(3); ...
                    ]));
                rep_volume.transforms = [rep_volume.transforms transform];
            end
        end
    end
    detector = petsird.DetectorModule(detecting_elements=rep_volume);
end

function crystal = get_crystal(cfg)
    % Return a cuboid shape
    crystal_shape = petsird.BoxShape(...
        corners = [...
            petsird.Coordinate(c=[0, 0, 0]), ...
            petsird.Coordinate(c=[0, 0, cfg.CRYSTAL_LENGTH(3)]), ...
            petsird.Coordinate(c=[0, cfg.CRYSTAL_LENGTH(2), cfg.CRYSTAL_LENGTH(3)]), ...
            petsird.Coordinate(c=[0, cfg.CRYSTAL_LENGTH(2), 0]), ...
            petsird.Coordinate(c=[cfg.CRYSTAL_LENGTH(1), 0, 0]), ...
            petsird.Coordinate(c=[cfg.CRYSTAL_LENGTH(1), 0, cfg.CRYSTAL_LENGTH(3)]), ...
            petsird.Coordinate(c=[cfg.CRYSTAL_LENGTH(1), cfg.CRYSTAL_LENGTH(2), cfg.CRYSTAL_LENGTH(3)]) ...
            petsird.Coordinate(c=[cfg.CRYSTAL_LENGTH(1), cfg.CRYSTAL_LENGTH(2), 0]), ...
        ]);
    crystal = petsird.BoxSolidVolume(shape=crystal_shape, material_id=1);
end

function efficiencies = get_detection_efficiencies(scanner, cfg)
    % Return some (non-physical) detection efficiencies
    % Only 1 type of module in the current scanner
    assert(length(scanner.scanner_geometry.replicated_modules) == 1);

    type_of_module = 1;
    num_detection_bins = petsird.helpers.get_num_detection_bins(scanner, type_of_module);
    event_energy_bin_edges = scanner.event_energy_bin_edges(type_of_module)(type_of_module);
    num_event_energy_bins = length(event_energy_bin_edges);

    detection_bin_efficiencies = ones(num_detection_bins, "single");

    rep_module = scanner.scanner_geometry.replicated_modules(1);
    num_modules = int32(length(rep_module.transforms));

    % We will use rotational symmetries translation along the axis
    % We assume all module-pairs are in coincidence, except those
    % with the same angle.
    % Writing a module number as (z-position, angle):
    %   eff((z1,a1), (z2, a2)) == eff((z1,0), (z2, abs(a2-a1)))
    % or in linear indices
    %   eff(z1 + NZ * a1, z2 + NZ * a2) == eff(z1, z2 + NZ * abs(a2 - a1))
    % (coincident) SGIDs need to start from 0, so ignoring self-coincident
    % angles we get
    %   SGID = z1 + NZ * (z2 + NZ * abs(a2 - a1) - 1)
    num_SGIDs = cfg.NUM_MODULES_ALONG_AXIS * cfg.NUM_MODULES_ALONG_AXIS * (cfg.NUM_MODULES_ALONG_RING - 1);
    NZ = cfg.NUM_MODULES_ALONG_AXIS;
    module_pair_SGID_LUT = zeros(num_modules, "int32");
    for mod1 = 0:num_modules-1
        for mod2 = 0:num_modules-1
            z1 = mod(mod1, NZ);
            a1 = idivide(mod1, NZ);
            z2 = mod(mod2, NZ);
            a2 = idivide(mod2, NZ);
            if a1 == a2
                module_pair_SGID_LUT(mod2+1, mod1+1) = -1;
            else
                module_pair_SGID_LUT(mod2+1, mod1+1) = z1 + NZ * (z2 + NZ * (abs(a2 - a1) - 1));
            end
        end
    end

    % fprint("SGID LUT:\n", module_pair_SGID_LUT, file=sys.stderr)
    assert(max(module_pair_SGID_LUT(:)) == num_SGIDs - 1);
    module_pair_efficiencies_vector = [];
    detecting_elements = rep_module.object.detecting_elements;
    num_det_els_in_module = length(detecting_elements.transforms);
    event_energy_bin_edges = scanner.event_energy_bin_edges(type_of_module);
    num_event_energy_bins = event_energy_bin_edges.number_of_bins();

    for SGID = 0:num_SGIDs-1
        % Extract first module_pair for this SGID. However, as this
        % currently unused, it is commented out
        % module_pair = numpy.argwhere(module_pair_SGID_LUT == SGID)[0]
        % print(module_pair, file=sys.stderr)
        module_pair_efficiencies = ones(...
                num_event_energy_bins * num_det_els_in_module, ...
                num_event_energy_bins * num_det_els_in_module ...
        );
        % give some (non-physical) value
        module_pair_efficiencies = module_pair_efficiencies * SGID;
        module_pair_efficiencies_vector = [module_pair_efficiencies_vector ...
            petsird.ModulePairEfficiencies(values=module_pair_efficiencies, sgid=SGID)];
        assert(length(module_pair_efficiencies_vector) == SGID + 1);
    end

    efficiencies = petsird.DetectionEfficiencies(...
        detection_bin_efficiencies=[detection_bin_efficiencies], ...
        module_pair_sgidlut=[[module_pair_SGID_LUT]], ...
        module_pair_efficiencies_vector=[[module_pair_efficiencies_vector]] ...
    );
end


function events = get_events(header, type_of_module_pair, num_events, cfg)
    % Generate some random events
    % TODO adapt code to use ExpandedDetectionBin as a better example for a user.
    % See the python/C++ code for how to do this.
    count1 = petsird.helpers.get_num_detection_bins(header.scanner, type_of_module_pair(1));
    count2 = petsird.helpers.get_num_detection_bins(header.scanner, type_of_module_pair(2));
    events = [];
    detection_bins = [petsird.DetectionBin(), petsird.DetectionBin()]
    for i = 1:num_events
        % Generate random detection_bins until detection effficiency is not zero
        while true
            detection_bins[1] = randi([0, count1-1]);
            detection_bins(2) = randi([0, count2-1]);
            if petsird.helpers.get_detection_efficiency(header.scanner, type_of_module_pair, detection_bins) > 0
                % in coincidence, we can get out of the loop
                break
            end
        end

        events = [events, petsird.CoincidenceEvent(...
                                                   detection_bins = detection_bins, ...
                                                   tof_idx=randi([0, cfg.NUMBER_OF_TOF_BINS-1])...
        )];
    end
end
