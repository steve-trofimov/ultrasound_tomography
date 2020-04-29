#include "argparse.h"

int main(int ac, char* av[]) {
    po::options_description desc("General options");
    std::string task_type;
    desc.add_options()
            ("help,h", "Show help")
            ("type,t", po::value<std::string>(&task_type), "Select task: ToF, detection, oneToF")
            ;
    po::options_description ToF_desc("ToF options");
    ToF_desc.add_options()
            ("input_dir,I", po::value<std::string>(), "Input directory with data")
            ("output,O", po::value<std::string>(), "Result ToF .txt file")
            ("info", po::value<bool>(), "Display information")
            ;
    po::options_description detection_desc("Detection options");
    detection_desc.add_options()
            ("input_water", po::value<std::string>(), "Input ToF for water")
            ("input_exp", po::value<std::string>(), "Input ToF for experiment")
            ("output,O", po::value<std::string>(), "Output matrix .txt file")
            ("info", po::value<bool>(), "Display information")
            ;
    po::options_description oneToF_desc("oneToF options");
    oneToF_desc.add_options()
            ("input_dir,I", po::value<std::string>(), "Input directory with data")
            ("number_emmiter", po::value<int>(), "Sensor number of interest")
            ("output_ToF,O", po::value<std::string>(), "Result ToF .txt file")
            ("output_read,O", po::value<std::string>(), "Read result .txt file")
            ;
    po::variables_map vm;
    try {
        po::parsed_options parsed = po::command_line_parser(ac, av).options(desc).allow_unregistered().run();
        po::store(parsed, vm);
        po::notify(vm);
        if(task_type == "ToF") {
            desc.add(ToF_desc);
            po::store(po::parse_command_line(ac, av, desc), vm);
            ToF(vm);
        }
        else if(task_type == "detection") {
            desc.add(detection_desc);
            po::store(po::parse_command_line(ac, av, desc), vm);
            detection(vm);
        }
        else if(task_type == "oneToF") {
            desc.add(oneToF_desc);
            po::store(po::parse_command_line(ac, av, desc), vm);
            oneToF(vm);
        }
        else {
            desc.add(ToF_desc).add(detection_desc);
            std::cout << desc << "\n";
            return -1;
        }
    }
    catch (std::exception& ex) {
        std::cout << desc << "\n";
        return -2;
    }
  return 0;
}