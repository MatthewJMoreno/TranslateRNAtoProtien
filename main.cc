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
#include "FastaParser.hh"
#include "SequenceTranslator.hh"
#include "FastaWriter.hh"
#include <iostream>
#include <cstdlib>
#include "timer.h"

int main(int argc, char *argv[])
{
	double readingTime;
	double computationTime;
	double writingTime;
	double totalTime;

	try
	{
		initialize_timer();

		const Parameters& parameters = get_parameters(argc, argv);

		std::auto_ptr<SequenceTranslator> sequence_translator = SequenceTranslatorFactory::create_translator(parameters.reading_frame, parameters.strand);
		FastaParser fasta_parser(parameters.input_file);
		FastaWriter fasta_writer(parameters.output_file);

		std::cout << "Reading FASTA nucleotide file..." << std::endl;
		start_timer();
		std::vector<FastaRecord> records = fasta_parser.parse_file();
		stop_timer();
   	readingTime = elapsed_time();
		std::cout << "Done. FASTA file has " << records.size() << " entries.\n" << std::endl;

		std::cout << "Translating to FASTA amino acid..." << std::endl;
		reset_timer();
		start_timer();
		sequence_translator->translate(records);
		stop_timer();
		computationTime = elapsed_time();
		std::cout << "Done\n" << std::endl;

		std::cout << "Writing FASTA amino acid file..." << std::endl;
		reset_timer();
		start_timer();
		fasta_writer.write(records);
		stop_timer();
		writingTime = elapsed_time();
		std::cout << "Done\n" << std::endl;

		totalTime = readingTime + computationTime + writingTime;

		std::cout << "Reading time: " << readingTime << 's' << '\n';
		std::cout << "Computation time: " << computationTime << 's' << '\n';
		std::cout << "Writing time: " << writingTime << 's' << '\n';
		std::cout << "Total time: " << totalTime << 's' << '\n';

		std::exit(0);
	}
	catch(const std::exception& e)
	{
		std::cerr << e.what() << std::endl;
		std::exit(-1);
	}
}
