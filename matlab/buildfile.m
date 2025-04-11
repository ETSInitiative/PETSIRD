function plan = buildfile
    plan = buildplan(localfunctions);
    plan.DefaultTasks = ["packageToolbox"];
end

function packageToolboxTask(~)
    buildToolbox("release");
end

function buildToolbox (outdir)
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end

    uuid = string(java.util.UUID.randomUUID());
    toolboxFolder = "./toolbox/";
    opts = matlab.addons.toolbox.ToolboxOptions(toolboxFolder, uuid);

    toolboxDirs = unique(fileparts(opts.ToolboxFiles));
    if ~all(contains(toolboxDirs, opts.ToolboxFolder))
        error("No symbolic links allowed in toolboxFolder (BUG through at least R2023a)");
    end

    opts.ToolboxName = "PETSIRD";

    opts.ToolboxVersion = "0.4.0";
    opts.OutputFile = fullfile(outdir, sprintf("petsird-%s.mltbx", opts.ToolboxVersion));

    opts.Description = "Positron Emission Tomography Standardization Initiative Raw Data (PETSIRD) toolbox for MATLAB";
    % opts.Summary = "";
    % opts.AuthorCompany = "";

    opts.MinimumMatlabRelease = "R2022b";

    % Must also specify which folders should be added to MATLAB path upon toolbox installation.
    % Must also include at least one *file *in the toolbox folder.
    % This seems to be bug on Linux, Matlab R2023a. On Windows, this isn't required.
    opts.ToolboxMatlabPath = toolboxFolder;

    matlab.addons.toolbox.packageToolbox(opts);
end
