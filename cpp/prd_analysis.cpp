#include "generated/hdf5/protocols.h"
#include <iostream>

int main(int argc, char *argv[])
{
    // Check if the user has provided a file
    if (argc < 2)
    {
        std::cerr << "Please provide a file to read" << std::endl;
        return 1;
    }

    // Open the file
    prd::hdf5::PrdExperimentReader reader(argv[1]);
    prd::Header header;
    reader.ReadHeader(header);

    std::cout << "Processing file: " << argv[1] << std::endl;
    std::cout << "Subject ID: " << header.subject.id << std::endl;
    std::cout << "Number of detectors: " << header.NumberOfDetectors() << std::endl;

    std::vector<prd::CoincidenceEvent> events;
    events.resize(100);

    // Process events in batches of up to 100
    float energy_1 = 0, energy_2 = 0;
    int num_events = 0;
    while (reader.ReadEvents(events))
    {
        for (auto &event : events)
        {
            energy_1 += static_cast<float>(event.energy_1);
            energy_2 += static_cast<float>(event.energy_2);
            num_events++;
        }
    }

    std::cout << "Number of events: " << num_events << std::endl;
    std::cout << "Average energy_1: " << energy_1 / num_events << std::endl;
    std::cout << "Average energy_2: " << energy_2 / num_events << std::endl;

    return 0;
}
