#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <sstream>
#include <cstdio>

void seq2tab(const std::string& input_table_file, const std::string& fasta_file, const std::string& id_column) {
    // Create a map to store sequences from the FASTA file
    std::unordered_map<std::string, std::string> sequences;

    // Read the FASTA file
    std::ifstream infile(fasta_file);
    if (!infile.is_open()) {
        std::cerr << "Error: Could not open FASTA file: " << fasta_file << std::endl;
        exit(1);
    }

    std::string line, current_id, current_sequence;
    while (std::getline(infile, line)) {
        if (line.length() > 0 && line[0] == '>') {
            if (!current_id.empty()) {
                sequences[current_id] = current_sequence;
            }
            // Remove '>' and extract ID before ';' if it exists
            std::string full_header = line.substr(1);
            size_t semicolon_pos = full_header.find(';');
            if (semicolon_pos != std::string::npos) {
                current_id = full_header.substr(0, semicolon_pos);
            } else {
                current_id = full_header;
            }
            current_sequence.clear();
        } else {
            current_sequence += line;
        }
    }
    if (!current_id.empty()) {
        sequences[current_id] = current_sequence;
    }
    infile.close();

    std::cout << "Loaded " << sequences.size() << " sequences from FASTA file." << std::endl;
    
    // Checkpoint: Print first sequence ID
    if (!sequences.empty()) {
        std::cout << "First sequence ID: " << sequences.begin()->first << std::endl;
    }

    // Read the input table and write the output table with the new sequence column
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

    std::string header;
    std::getline(table_infile, header);
    
    // Checkpoint: Check if id_column exists in header
    std::istringstream header_stream(header);
    std::string column_name;
    int id_column_index = -1;
    int col_index = 0;
    while (std::getline(header_stream, column_name, '\t')) {
        if (column_name == id_column) {
            id_column_index = col_index;
            std::cout << "Found id_column '" << id_column << "' at position " << col_index << std::endl;
            break;
        }
        col_index++;
    }
    
    if (id_column_index == -1) {
        std::cerr << "Warning: Column '" << id_column << "' not found in header. Using first column instead." << std::endl;
        id_column_index = 0;
    }
    
    table_outfile << header << "\tsequence\n"; // Use tab as the delimiter

    int matched = 0, total = 0;
    std::string row;
    bool first_row = true;
    while (std::getline(table_infile, row)) {
        std::istringstream row_stream(row);
        std::string cell, id;
        
        // Extract the ID from the specified column
        for (int i = 0; i <= id_column_index; i++) {
            std::getline(row_stream, id, '\t');
        }
        
        // Checkpoint: Print first value from id_column
        if (first_row) {
            std::cout << "First value in id_column: " << id << std::endl;
            first_row = false;
        }

        std::string sequence = "NA";
        if (sequences.find(id) != sequences.end()) {
            sequence = sequences[id];
            matched++;
        }
        
        table_outfile << row << "\t" << sequence << "\n"; // Use tab as the delimiter
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

    std::cout << "Successfully added sequences to table." << std::endl;
    std::cout << "Matched " << matched << " out of " << total << " rows." << std::endl;
}

int main(int argc, char* argv[]) {
    // Check command line arguments
    if (argc != 4) {
        std::cout << "Usage: " << argv[0] << " <input_table_file> <fasta_file> <id_column>" << std::endl;
        std::cout << "Example: " << argv[0] << " data.tsv sequences.fasta id" << std::endl;
        std::cout << std::endl;
        std::cout << "This program adds a 'sequence' column to a tab-separated table file" << std::endl;
        std::cout << "by matching IDs from the first column with sequence headers in a FASTA file." << std::endl;
        std::cout << "The original table file will be modified in-place." << std::endl;
        return 1;
    }

    std::string input_table_file = argv[1];
    std::string fasta_file = argv[2];
    std::string id_column = argv[3]; // Currently not used, assumes first column

    // Check if input files exist
    std::ifstream test_table(input_table_file);
    if (!test_table.good()) {
        std::cerr << "Error: Input table file '" << input_table_file << "' does not exist or cannot be read." << std::endl;
        return 1;
    }
    test_table.close();

    std::ifstream test_fasta(fasta_file);
    if (!test_fasta.good()) {
        std::cerr << "Error: FASTA file '" << fasta_file << "' does not exist or cannot be read." << std::endl;
        return 1;
    }
    test_fasta.close();

    // Call the seq2tab function
    seq2tab(input_table_file, fasta_file, id_column);
    
    return 0;
}