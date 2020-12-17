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

#include <iostream>
#include "FastaParser.hh"
#include "FastaWriter.hh"
#include <cstdlib>
#include "timer.h"
#include <filesystem>
#include "mpi.h"
#include "Parameters.hh"
#include "SequenceTranslator.hh"
#include <fstream>
#include "FastaRecord.hh"

void MPI_Parallelization(Parameters parameters, char* fileName, int id, int p, int fileWidth, long long fileSize){
	std::vector<FastaRecord> records;
	std::auto_ptr<SequenceTranslator> sequence_translator = SequenceTranslatorFactory::create_translator(parameters.reading_frame, parameters.strand);

	MPI_Offset globalstart;
	MPI_File in, out;
	char* nucleotides;
	char* headers = (char*)malloc(sizeof(p));
	long bufferSize, offset, actualEnd, writeOffset = 0;
	long *allWriteOffsets = (long*)malloc(sizeof(long) * p);
	long actualStart = 0;
	double localReadinTime, globalReadinTime, localComputationTime, globalComputationTime, localWriteInTime, globalWriteInTime, totalTime;
	initialize_timer();

	//Each processor has a buffer of size offset + filewidth
	offset = fileSize/p;
	bufferSize = offset + fileWidth + 1;
	nucleotides = (char*)malloc(bufferSize);

	if (nucleotides == NULL){
		std::cout << "failed to allocate memory" << std::endl;
	}

	//Determine where each node starts their read
	globalstart = offset * id;

	if (p == 1){
		actualEnd = fileSize;
		//std::cout << actualEnd << std::endl;
	} else {
		actualEnd = offset;
	}

	offset += fileWidth;

	//Now read in the file
	start_timer();
	std::cout << "Process " << id << " is reading in the file" << std::endl;
	int failedToOpen = MPI_File_open(MPI_COMM_WORLD, fileName, MPI_MODE_RDONLY, MPI_INFO_NULL, &in);
	if (failedToOpen){
		std::cout << "Failed to open: " << fileName << std::endl;
	}

	if (p == 1){
		failedToOpen = MPI_File_read(in, nucleotides, fileSize, MPI_CHAR, MPI_STATUS_IGNORE);
	} else {
		failedToOpen = MPI_File_read_at(in, globalstart, nucleotides, offset, MPI_CHAR, MPI_STATUS_IGNORE);
	}

	if (failedToOpen){
		std::cout << "Can't read from " << fileName << " this may be because its too big try using more machines." << std::endl;
		std::exit(0);
	}

	//Determine the actual starting and ending positions each node is responsible for.
	//Some nodes have more housekeeping to do than others
	//node 0 only needs to make sure its actual end is correct
	if (id == 0 && p != 1){
		while(true){
			if (nucleotides[actualEnd] != '\n'){
				actualEnd++;
			} else {
				break;
			}
		}
	}

	//Middle nodes needs to update their actualStart and actualEnd
	if (id > 0 && id < p - 1){
		//update start
		while(true){
			if (nucleotides[actualStart] != '\n'){
				actualStart++;
			} else {
				break;
			}
		}

		//update end
		while(true){
			if (nucleotides[actualEnd] != '\n'){
				actualEnd++;
			} else {
				break;
			}
		}
	}

	//The last node only needs to update its start
	if (id == p-1 && p != 1){
		//update start
		while(true){
			if (nucleotides[actualStart] != '\n'){
				actualStart++;
			} else {
				break;
			}
		}
	}

	//Now we extract the data out of the char buffer nucleotides and create FastaRecords to add to records
	//The computation requires a vector of FastaRecord
	while(actualStart < actualEnd){
		FastaRecord newRecord = FastaRecord();
		//>$ serves as an indicator that this record is part of a different nodes data
		newRecord.header = ">$";
		newRecord.sequence.resize(offset);
		long sequenceIndex = 0;
		int stillReadingRecord = 1;

		while (stillReadingRecord == 1){
			if (sequenceIndex == 0 && (nucleotides[actualStart] == '\n')){
				actualStart++;
				continue;
			} else {
				if (nucleotides[actualStart] == '>'){
					int index = 0;
					newRecord.header = "";
					newRecord.header.resize(20);
					while (nucleotides[actualStart] != '\n'){
						newRecord.header[index] = nucleotides[actualStart];
						index++;
						actualStart++;
					}
					newRecord.header.resize(index);
				} else {
					while (true){
						if (nucleotides[actualStart] == '>' || actualStart >= actualEnd){
							stillReadingRecord = 0;
							break;
						} else {
							if (nucleotides[actualStart] == '\n'){
								actualStart++;
							} else {
								newRecord.sequence[sequenceIndex] = nucleotides[actualStart];
								sequenceIndex++;
								actualStart++;
							}
						}
					}
				}
			}
		}

		newRecord.sequence.resize(sequenceIndex);
		records.push_back(newRecord);
	}

	stop_timer();
	localReadinTime = elapsed_time();
	std::cout << "Processor " << id << " is done with reading in the data and took " << localReadinTime << " seconds" << std::endl;

	std::cout << "Processor " << id << " has started the computation " << std::endl;
	reset_timer();
	start_timer();
	sequence_translator->translate(records);
	stop_timer();
	localComputationTime = elapsed_time();
	std::cout << "Processor " << id << " is done with the computation and took " << localComputationTime << " seconds" << std::endl;

	std::cout << "Processor " << id << " has started writing to the file" << std::endl;
	reset_timer();
	start_timer();
	long totalChars = 0;
	long index = 0;
	std::string data;
	data.resize(offset);
	for (unsigned long i = 0; i < records.size(); i++){
		//Only copy the header in if its not equal to partial
		if (records[i].header.compare(">$") != 0){
			for (unsigned long j = 0; j < records[i].header.length(); j++){
				data[index] = records[i].header[j];
				totalChars++;
				index++;
			}

			data[index] = '\n';
			totalChars++;
			index++;
		}

		//Copy the sequence into the data string
		for (unsigned long j = 0; j < records[i].sequence.length(); j++){
			data[index] = records[i].sequence[j];
			index++;
			totalChars++;
		}

		//Check if we are done with the partial record, if so add a newline char at the end
		if (i+1 < records.size() && records[i+1].header.compare(">$") != 0){
			data[index] = '\n';
			totalChars++;
			index++;
		}
	}

	//Determine if process i+1's first header is >$ if so we don't add a newline character to the end of the last sequence in records
	char firstHeader = records[0].header[1];
	MPI_Allgather(&firstHeader, 1, MPI_CHAR, headers, 1, MPI_CHAR, MPI_COMM_WORLD);
	if (id < p - 1){
		if (headers[id + 1] != '$'){
			data[index] = '\n';
			totalChars++;
		}
	}

	data.resize(totalChars);

	long stringLen = (long)data.length();
	MPI_Allgather(&stringLen, 1, MPI_LONG, allWriteOffsets, 1, MPI_LONG, MPI_COMM_WORLD);

	//The actual writing offset of process i is determined by the stringLen of process i-1 + stringlen of process i-2 + ... + stringlen of process 0
	for (int i = 0; i < id; i++){
		writeOffset+=allWriteOffsets[i];
	}

	MPI_File_open(MPI_COMM_WORLD, parameters.output_file.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &out);
	MPI_File_write_at(out, writeOffset, data.c_str(), stringLen, MPI_CHAR, MPI_STATUS_IGNORE);
	stop_timer();
	localWriteInTime = elapsed_time();

	std::cout << "Processor " << id << " has finished writing to the file and took " << localWriteInTime << std::endl;

	free(headers);
	free(nucleotides);
	MPI_File_close(&in);
	MPI_File_close(&out);

	MPI_Reduce(&localReadinTime, &globalReadinTime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&localComputationTime, &globalComputationTime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&localWriteInTime, &globalWriteInTime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (id == 0){
		totalTime = globalWriteInTime + globalReadinTime + globalComputationTime;
		std::cout << "Final reading time: " << globalReadinTime / p << " seconds" << std::endl;
		std::cout << "Final computation time: " << globalComputationTime / p << " seconds" << std::endl;
		std::cout << "Final writing time: " << globalWriteInTime / p << " seconds" << std::endl;
		std::cout << "Final total time: " << totalTime / p << " seconds" << std::endl;
	}
}

int main(int argc, char *argv[]){
	int id;
	int p;
	long fileWidth;
	long fileSize;
	std::string line;

	try{
		const Parameters& parameters = get_parameters(argc, argv);
		FastaWriter fasta_writer(parameters.output_file);

		MPI_Init( &argc, &argv );
		MPI_Comm_rank( MPI_COMM_WORLD, &id ); //rank is the processes position in the order
		MPI_Comm_size( MPI_COMM_WORLD, &p ); //Determine the total number of processes in communicator

		//Find the first sequence in the file, this will be used to set the width of the file
		std::ifstream ifs;
		ifs.open(argv[1]);
		while (std::getline(ifs, line)) {
			if (line[0] !=  '>') {
				fileWidth = line.size();
				break;
			}
		}
		ifs.close();

		std::filesystem::path file{argv[1]};
		fileSize = std::filesystem::file_size(file);

		MPI_Parallelization(parameters, argv[1], id, p, fileWidth, fileSize);
	} catch(const std::exception& e) {
		std::cerr << e.what() << std::endl;
	}

	MPI_Finalize();
	return 0;
}
