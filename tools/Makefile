# Makefile for rename_fasta standalone executable

CXX = g++
CXXFLAGS = -std=c++11 -O2 -Wall

PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin

# Targets
RENAME_TARGET = rename_fasta
RENAME_SOURCE = rename_fasta_standalone.cpp

SEQ2TAB_TARGET = seq2tab
SEQ2TAB_SOURCE = seq2tab_standalone.cpp

ASSIGN_TARGET = assign_motu
ASSIGN_SOURCE = assign_motu_standalone.cpp

# Default target - build both
all: $(RENAME_TARGET) $(SEQ2TAB_TARGET) $(ASSIGN_TARGET)

# Individual targets
$(RENAME_TARGET): $(RENAME_SOURCE)
	$(CXX) $(CXXFLAGS) -o $(RENAME_TARGET) $(RENAME_SOURCE)

$(SEQ2TAB_TARGET): $(SEQ2TAB_SOURCE)
	$(CXX) $(CXXFLAGS) -o $(SEQ2TAB_TARGET) $(SEQ2TAB_SOURCE)

$(ASSIGN_TARGET): $(ASSIGN_SOURCE)
	$(CXX) $(CXXFLAGS) -o $(ASSIGN_TARGET) $(ASSIGN_SOURCE)

# Clean target
clean:
	rm -f $(RENAME_TARGET) $(SEQ2TAB_TARGET) $(ASSIGN_TARGET)
# Install target (optional - copies to /usr/local/bin)
install: $(RENAME_TARGET) $(SEQ2TAB_TARGET) $(ASSIGN_TARGET)
	cp $(RENAME_TARGET) $(BINDIR)/$(RENAME_TARGET)
	cp $(SEQ2TAB_TARGET) $(BINDIR)/$(SEQ2TAB_TARGET)
	cp $(ASSIGN_TARGET) $(BINDIR)/$(ASSIGN_TARGET)
# Uninstall target (optional)
uninstall:
	rm -f $(BINDIR)/$(RENAME_TARGET) $(BINDIR)/$(SEQ2TAB_TARGET) $(BINDIR)/$(ASSIGN_TARGET)
.PHONY: clean install uninstall