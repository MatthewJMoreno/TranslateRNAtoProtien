/*
 * (C) Copyright 2014-2015 Enterome
 * 
 * This file is part of sequence-translator.
 * 
 * sequence-translator is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * sequence-translator is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with sequence-translator.  If not, see <http://www.gnu.org/licenses/>.
 */ 

#include "Parameters.hh"
#include <iostream>
#include <fstream>
#include <stdexcept>

void help_message(){
	std::cout << "\nThis program takes the following arguments:" << '\n';
	std::cout << "  inputFile: A FASTA nucleotide file to convert" << '\n';
	std::cout << "  outputFile: FASTA amino acid file generated" << '\n';
	std::cout << "  reading-frame *optional*: 1, 2 or 3. Note, default is set to 1" << '\n';
	std::cout << "  strand *optional*: direct or reverse. Note, default is set to direct" << '\n' << '\n';
}

Parameters get_parameters(int argc, char* argv[])
{
	Parameters pars;

	if (argc < 3){
		help_message();
		std::exit(0);
	}

	pars.input_file = argv[1];
	pars.output_file = argv[2];
	pars.reading_frame = 1;
	pars.strand = "direct";

	if (argc > 3){
		pars.reading_frame = atoi(argv[3]);
		if (pars.reading_frame < 1 || pars.reading_frame > 3){
			help_message();
			std::exit(0);
		}
	}

	if (argc > 4){
		pars.strand = argv[4];
		if (!(pars.strand.compare("direct") == 0 || pars.strand.compare("reverse") == 0)){
			help_message();
			std::exit(0);
		}
	}

	if (argc > 5){
		help_message();
		std::exit(0);
	}

	std::ifstream ifs;
	ifs.open(pars.input_file.c_str());

	if (ifs.fail())
	{
		throw (std::invalid_argument("error: " + pars.input_file + 
					" cannot be opened. Check that the path is valid and that you have read permissions."));
	}

	ifs.close();


	// Check output file
	std::ofstream ofs;
	ofs.open(pars.output_file.c_str());

	if (ofs.fail())
	{
		throw (std::invalid_argument("error: " + pars.output_file + 
					" cannot be created. Check that the path is valid and that you have write permissions."));
	}

	ofs.close();

	return (pars);
}
