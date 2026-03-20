#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <sstream>
#include <cstdio>
#include <vector>

void assign_motu(const std::string& swarm_file, const std::string& input_table_file) {
    // Create a map to store MOTU assignments from SWARM output
    // Key: sequence ID, Value: MOTU name (representative sequence)
    std::unordered_map<std::string, std::string> seq_to_motu;

    // Read the SWARM output file
    std::ifstream swarm_infile(swarm_file);
    if (!swarm_infile.is_open()) {
        std::cerr << "Error: Could not open SWARM output file: " << swarm_file << std::endl;
        exit(1);
    }

    std::string line;
    int total_motus = 0;
    int total_sequences = 0;
    
    while (std::getline(swarm_infile, line)) {
        if (line.empty()) continue;
        
        total_motus++;
        
        // Split line by space to get individual sequence headers
        std::vector<std::string> sequence_headers;
        std::istringstream line_stream(line);
        std::string seq_header;
        
        while (line_stream >> seq_header) {
            if (!seq_header.empty()) {
                sequence_headers.push_back(seq_header);
            }
        }
        
        if (sequence_headers.empty()) continue;
        
        // The first sequence is the MOTU representative (seed)
        // Extract the ID from the first header (before semicolon if it exists)
        std::string motu_name = sequence_headers[0];
        size_t semi_pos = motu_name.find(';');
        if (semi_pos != std::string::npos) {
            motu_name = motu_name.substr(0, semi_pos);
        }
        
        // Assign all sequences in this cluster to this MOTU
        for (const std::string& seq_header : sequence_headers) {
            // Extract just the sequence ID (before any semicolon)
            std::string seq_id = seq_header;
            size_t semi_pos = seq_id.find(';');
            if (semi_pos != std::string::npos) {
                seq_id = seq_id.substr(0, semi_pos);
            }
            
            seq_to_motu[seq_id] = motu_name;
            total_sequences++;
        }
    }
    swarm_infile.close();

    std::cout << "Loaded " << total_motus << " MOTUs from SWARM output." << std::endl;
    std::cout << "Total sequences mapped: " << total_sequences << std::endl;

    // Read the input table and write the output table with the MOTU column
    std::ifstream table_infile(input_table_file);
    std::string temp_table_file = input_table_file + ".tmp";
    std::ofstream table_outfile(temp_table_file);
    
    if (!table_infile.is_open()) {
        std::cerr << "Error: Could not open input table file: " << input_table_file << std::endl;
        exit(1);
    }
    
    if (!table_outfile.is_open()) {
        std::cerr << "Error: Could not open temporary output table file: " << temp_table_file << std::endl;
        exit(1);
    }

    // Read and modify header
    std::string header;
    std::getline(table_infile, header);
    
    // Find the position to insert MOTU column (after ID column)
    std::istringstream header_stream(header);
    std::vector<std::string> header_cols;
    std::string col;
    
    while (std::getline(header_stream, col, '\t')) {
        header_cols.push_back(col);
    }
    
    // Write header with MOTU column after ID
    table_outfile << header_cols[0] << "\tMOTU";
    for (size_t i = 1; i < header_cols.size(); i++) {
        table_outfile << "\t" << header_cols[i];
    }
    table_outfile << "\n";

    int matched = 0, unmatched = 0, total = 0;
    std::string row;
    
    while (std::getline(table_infile, row)) {
        if (row.empty()) continue;
        
        std::istringstream row_stream(row);
        std::string id;
        std::getline(row_stream, id, '\t');
        
        // Get remaining columns
        std::string remaining_cols;
        std::getline(row_stream, remaining_cols);
        
        // Find MOTU assignment
        std::string motu = "NA";
        if (seq_to_motu.find(id) != seq_to_motu.end()) {
            motu = seq_to_motu[id];
            matched++;
        } else {
            unmatched++;
        }
        
        // Write row with MOTU column
        table_outfile << id << "\t" << motu << "\t" << remaining_cols << "\n";
        total++;
    }

    table_infile.close();
    table_outfile.close();

    // Replace the original file with the temporary file
    if (std::remove(input_table_file.c_str()) != 0) {
        std::cerr << "Error: Could not remove original file: " << input_table_file << std::endl;
        exit(1);
    }
    
    if (std::rename(temp_table_file.c_str(), input_table_file.c_str()) != 0) {
        std::cerr << "Error: Could not rename temporary file to: " << input_table_file << std::endl;
        exit(1);
    }

    std::cout << "Successfully assigned MOTUs to table." << std::endl;
    std::cout << "Matched " << matched << " out of " << total << " rows." << std::endl;
    if (unmatched > 0) {
        std::cout << "Warning: " << unmatched << " sequences could not be assigned to any MOTU." << std::endl;
    }
}

int main(int argc, char* argv[]) {
    // Check command line arguments
    if (argc != 3) {
        std::cout << "Usage: " << argv[0] << " <swarm_output> <input_table>" << std::endl;
        std::cout << "Example: " << argv[0] << " EXPX_ODIN_SWARM_output EXPX_ODIN_seqs.csv" << std::endl;
        std::cout << std::endl;
        std::cout << "This program assigns MOTU names to sequences in a table based on SWARM clustering output." << std::endl;
        std::cout << "The SWARM output file should contain clusters with sequences separated by '; '." << std::endl;
        std::cout << "The first sequence in each cluster is used as the MOTU representative name." << std::endl;
        std::cout << "A new 'MOTU' column will be added after the ID column." << std::endl;
        std::cout << "The original table file will be modified in-place." << std::endl;
        return 1;
    }

    std::string swarm_file = argv[1];
    std::string input_table_file = argv[2];

    // Check if input files exist
    std::ifstream test_swarm(swarm_file);
    if (!test_swarm.good()) {
        std::cerr << "Error: SWARM output file '" << swarm_file << "' does not exist or cannot be read." << std::endl;
        return 1;
    }
    test_swarm.close();

    std::ifstream test_table(input_table_file);
    if (!test_table.good()) {
        std::cerr << "Error: Input table file '" << input_table_file << "' does not exist or cannot be read." << std::endl;
        return 1;
    }
    test_table.close();

    // Call the assign_motu function
    assign_motu(swarm_file, input_table_file);
    
    return 0;
}
