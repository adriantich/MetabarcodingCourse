#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>

void rename_fasta(const std::string& input_file, const std::string& output_file, const std::string& prefix) {
    std::ifstream infile(input_file);
    std::ofstream outfile(output_file);
    std::string line;
    int count = 1;

    if (!infile.is_open()) {
        std::cerr << "Error: Could not open input file: " << input_file << std::endl;
        exit(1);
    }
    
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open output file: " << output_file << std::endl;
        exit(1);
    }

    while (std::getline(infile, line)) {
        if (line.length() > 0 && line[0] == '>') {
            std::ostringstream new_header;
            new_header << ">" << prefix << "_" << std::setw(9) << std::setfill('0') << count;
            
            // Find the first semicolon and preserve everything after it
            size_t semicolon_pos = line.find(';');
            if (semicolon_pos != std::string::npos) {
                new_header << line.substr(semicolon_pos);
            }
            
            outfile << new_header.str() << std::endl;
            count++;
        } else {
            outfile << line << std::endl;
        }
    }

    infile.close();
    outfile.close();
    
    std::cout << "Successfully renamed FASTA sequences in '" << input_file 
              << "' to '" << output_file << "' with prefix '" << prefix << "'" << std::endl;
}

int main(int argc, char* argv[]) {
    // Check command line arguments
    if (argc != 4) {
        std::cout << "Usage: " << argv[0] << " <input_file> <output_file> <prefix>" << std::endl;
        std::cout << "Example: " << argv[0] << " input.fasta output.fasta SEQ" << std::endl;
        return 1;
    }

    std::string input_file = argv[1];
    std::string output_file = argv[2];
    std::string prefix = argv[3];

    // Check if input file exists
    std::ifstream test_file(input_file);
    if (!test_file.good()) {
        std::cerr << "Error: Input file '" << input_file << "' does not exist or cannot be read." << std::endl;
        return 1;
    }
    test_file.close();

    // Call the rename function
    rename_fasta(input_file, output_file, prefix);
    
    return 0;
}