#include <iostream>

#include "generated/hdf5/protocols.h"

// return pair of integers between 0 and max
std::pair<int, int> get_random_pair(int max)
{
    int a = rand() % max;
    int b = rand() % max;
    return std::make_pair(a, b);
}

// Random float between 0 and 1
float get_random_float()
{
    return (float)rand() / (float)RAND_MAX;
}

int main()
{
    prd::Header head;

    // Create some detectors
    std::vector<float> offsets = {0.0, 1.0, 2.0};
    std::vector<float> angles;
    for (int i = 0; i < 10; ++i)
    {
        angles.push_back(2 * M_PI * i / 10);
    }

    int detector_id = 0;
    for (auto offset : offsets)
    {
        for (auto angle : angles)
        {
            // Create a new detector
            prd::Detector d;
            d.offset = offset;
            d.angle = angle;
            d.id = detector_id++;
            head.detectors.push_back(d);
        }
    }

    head.institution.name = "Diamond Light Source";
    head.institution.address = "Harwell Science and Innovation Campus, Didcot, Oxfordshire, OX11 0DE, UK";
    head.subject.id = "123456";

    std::string outfile = "test.h5";
    prd::hdf5::PrdExperimentWriter writer(outfile);
    writer.WriteHeader(head);

    // Create some events
    const int num_events = 1000;
    std::vector<prd::CoincidenceEvent> events;
    for (int i = 0; i < num_events; ++i)
    {
        auto detectors = get_random_pair(head.NumberOfDetectors());
        prd::CoincidenceEvent e;
        e.detector1_id = detectors.first;
        e.detector2_id = detectors.second;
        e.energy1 = get_random_float();
        e.energy2 = get_random_float();
        e.delta_t = get_random_float();
        events.push_back(e);
    }
    // Illustrates writing events from a buffer rather than individual events
    writer.WriteSamples(events);
    writer.EndSamples();

    // Check that we have completed protocol
    writer.Close();
    return 0;
}
