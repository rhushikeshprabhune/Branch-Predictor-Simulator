#include <iostream>
#include <string.h>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <iomanip>
using namespace std;

typedef struct{
unsigned long t_tag, t_LRU_counter, t_dirty_bit, t_valid_bit; //Define structure for BTB
}btb_table; 

int main(int argc, char * argv[])
{
	char * trace_file;
	unsigned long M1, M2, N, K, BTB_SIZE, BTB_ASSOC;
	char * predictor = argv[1];
	//cmd line args for bimodal branch predictor
	if (!strcmp(argv[1], "bimodal"))
	{
		if (argc == 6)
		{
			M2 = atoi(argv[2]);
			BTB_SIZE = atoi(argv[3]);
			BTB_ASSOC = atoi(argv[4]);
			trace_file = argv[5];
		}
		else
		{
			cout << "Wrong number of arguments!" << endl;
		}
	}
	else if (!strcmp(argv[1], "gshare")) //cmd line args for gshare branch predictor
	{
		if (argc == 7)
		{
			M1 = atoi(argv[2]);
			N = atoi(argv[3]);
			BTB_SIZE = atoi(argv[4]);
			BTB_ASSOC = atoi(argv[5]);
			trace_file = argv[6];
		}
		else
		{
			cout << "Wrong number of arguments!" << endl;
		}
	}
	else if (!strcmp(argv[1], "hybrid")) //cmd line args for hybrid branch predictor
	{
		if (argc == 9)
		{
			K = atoi(argv[2]);
			M1 = atoi(argv[3]);
			N = atoi(argv[4]);
			M2 = atoi(argv[5]);
			BTB_SIZE = atoi(argv[6]);
			BTB_ASSOC = atoi(argv[7]);
			trace_file = argv[8];
		}
		else
		{
			cout << "Wrong number of arguments!" << endl;
		}

	}


	unsigned long BLOCKSIZE, sets, block, bo, indexsize, tagsize;
	unsigned long reads, writes, read_misses, write_misses;
	BLOCKSIZE = 4;
	block = BTB_SIZE / BLOCKSIZE; //number of blocks
	if (BTB_ASSOC != 0 && BTB_SIZE != 0)
	{
		sets = BTB_SIZE / (BLOCKSIZE * BTB_ASSOC);
	}
	else
	{
		sets = 1; //number of sets
	}
	bo = log2(BLOCKSIZE); //block offset
	indexsize = log2(sets); //number of bits in index value
	tagsize = 24 - bo - indexsize; //tag length                                
	reads = 0;
	writes = 0;
	read_misses = 0;
	write_misses = 0;
	btb_table ** btb_cache;
	btb_cache = new btb_table * [sets]; //Define the two dimentional array for BTB table
	for (unsigned int i = 0; i < sets; i++)
	{
		btb_cache[i] = new btb_table[BTB_ASSOC];
	}

	for (unsigned int i=0;i<sets;i++)
	{
		for (unsigned int j = 0; j < BTB_ASSOC; j++)
		{
			btb_cache[i][j].t_LRU_counter = j; //initialise the lru counters, dirty bits,valid bits and tags
			btb_cache[i][j].t_dirty_bit = 0;
			btb_cache[i][j].t_valid_bit = 0;
			btb_cache[i][j].t_tag = 0;
		}
	}


	unsigned long bimodal_table_size = pow(2, M2); //declare prediction table for bimodal
	unsigned int * pred_table_bimodal;
	pred_table_bimodal = new unsigned int[bimodal_table_size];
	for (unsigned long i = 0; i < bimodal_table_size; i++)
	{
		pred_table_bimodal[i] = 2;
	}

	unsigned long gshare_table_size = pow(2, M1); //declare prediction table for gshare
	unsigned int * pred_table_gshare;
	pred_table_gshare = new unsigned int[gshare_table_size];
	for (unsigned int j = 0; j < gshare_table_size; j++)
	{
		pred_table_gshare[j] = 2;
	}

	unsigned long chooser_table_size = pow(2, K); //chooser table for hybrid
	unsigned int * chooser_table;
	chooser_table = new unsigned int[chooser_table_size];
	for (unsigned int j = 0; j < chooser_table_size; j++)
	{
		chooser_table[j] = 1;
	}

	unsigned long misprediction = 0;
	unsigned long num_predictions = 0;
	unsigned long prediction = 0;
	unsigned long reality = 0;
	unsigned long bhr = 0;
	unsigned long prediction_bimodal = 0;
	unsigned long prediction_gshare = 0;
	unsigned long proceedtoBP = 0;
	unsigned long num_mispred_bp = 0;
	unsigned long num_pred_bp = 0;
	unsigned long num_btb_missandtaken = 0;
	unsigned long btb_miss = 0;
	int isitbimodal = 0;
	int isitgshare = 0;
	int isithybrid = 0;

	ifstream infile(trace_file);
	unsigned long branch;
	char takenornot;
	while (infile >> hex >> branch >> takenornot) //read the branch address and branch result from trace file
	{
		num_predictions++;
		if (!strcmp(predictor, "bimodal"))
		{
			isitbimodal = 1;
			unsigned long index_bimodal = ((1 << M2) - 1) & (branch >> 2);
			if (BTB_SIZE != 0 && BTB_ASSOC != 0)
			{
				unsigned long tag = branch >> (bo + indexsize); //extract the tag from input address
				unsigned long btb_index = ((1 << indexsize) - 1) & (branch >> (bo));
				int r_hit_block = 0;
				int BTB_hit = 0;
				for (unsigned long j = 0; j < BTB_ASSOC; j++) //find if hit or miss and if hit find the hit block   
				{
					if (btb_cache[btb_index][j].t_tag == tag)
					{
						BTB_hit = 1;
						r_hit_block = j;
					}
				}
				if (BTB_hit == 1)
				{ //hit case
					unsigned int oldLRU;
					oldLRU = btb_cache[btb_index][r_hit_block].t_LRU_counter;
					for (unsigned long i = 0; i < BTB_ASSOC; i++)
					{
						if (btb_cache[btb_index][i].t_LRU_counter < oldLRU) //update lru counters if the value of counters
							btb_cache[btb_index][i].t_LRU_counter += 1; //is less than hit block counter value
					}
					btb_cache[btb_index][r_hit_block].t_LRU_counter = 0; //put lru counter of hit block to zero
					proceedtoBP = 1;
				}


				if (BTB_hit == 0)
				{ //miss case 
					read_misses++;
					int invalid = 0;
					prediction = 0;
					btb_miss = 1;

					for (unsigned int vl = 0; vl < BTB_ASSOC; vl++)
					{
						if (btb_cache[btb_index][vl].t_valid_bit == 0)
						{ //check for number of invalids
							invalid++;
						}
					}


					if (invalid == 1)
					{ //case when only one block is empty
						for (unsigned int i = 0; i < BTB_ASSOC; i++)
						{
							if (btb_cache[btb_index][i].t_valid_bit == 0)
							{ //update tag,LRU_counter and valid_bit 
								for (unsigned int j = 0; j < BTB_ASSOC; j++)
								{
									btb_cache[btb_index][j].t_LRU_counter++;
								}

								btb_cache[btb_index][i].t_tag = tag;
								btb_cache[btb_index][i].t_LRU_counter = 0;
								btb_cache[btb_index][i].t_valid_bit = 1;
							}
						}
					}

					if (invalid == 0)
					{ //when all are filled, check highest LRU_count and update tag to that block
						unsigned int temp, highest_LRU;
						temp = btb_cache[btb_index][0].t_LRU_counter;
						highest_LRU = 0;
						for (unsigned int j = 0; j < BTB_ASSOC; j++)
						{
							if (btb_cache[btb_index][j].t_LRU_counter > temp)
							{
								temp = btb_cache[btb_index][j].t_LRU_counter;
								highest_LRU = j;
							}
						}

						for (unsigned int k = 0; k < BTB_ASSOC; k++)
						{
							btb_cache[btb_index][k].t_LRU_counter++;
						}

						btb_cache[btb_index][highest_LRU].t_tag = tag;
						btb_cache[btb_index][highest_LRU].t_LRU_counter = 0;
						btb_cache[btb_index][highest_LRU].t_dirty_bit = 0;

					}

					if (invalid != 0 && invalid != 1)
					{ //case when multiple blocks are empty     
						unsigned int currentHighest = 0, highest_LRU;
						for (unsigned int k = 0; k < BTB_ASSOC; k++)
						{
							if (btb_cache[btb_index][k].t_valid_bit == 0 && btb_cache[btb_index][k].t_LRU_counter > currentHighest)
							{ //check which are empty and find the one with highest lru count
								currentHighest = btb_cache[btb_index][k].t_LRU_counter;
								highest_LRU = k;
							}

						}

						//Update the LRU Counters
						for (unsigned int k = 0; k < BTB_ASSOC; k++)
						{
							if (btb_cache[btb_index][k].t_LRU_counter < btb_cache[btb_index][highest_LRU].t_LRU_counter)
								btb_cache[btb_index][k].t_LRU_counter++;
						}

						btb_cache[btb_index][highest_LRU].t_tag = tag; //update tag,lru counter and valid bit 
						btb_cache[btb_index][highest_LRU].t_LRU_counter = 0;
						btb_cache[btb_index][highest_LRU].t_valid_bit = 1;
					}

				}
			}
			//start bimodal logic
			if (proceedtoBP == 1 || (BTB_ASSOC == 0 && BTB_SIZE == 0))
			{
				num_pred_bp++;
				if ((pred_table_bimodal[index_bimodal] >> 1) == 1)
					prediction = 1;
				else if ((pred_table_bimodal[index_bimodal] >> 1) == 0)
					prediction = 0;

				if (takenornot == 't')
				{
					if (pred_table_bimodal[index_bimodal] != 3)
						pred_table_bimodal[index_bimodal] = pred_table_bimodal[index_bimodal] + 1;
				}
				else if (takenornot == 'n')
				{
					if (pred_table_bimodal[index_bimodal] != 0)
						pred_table_bimodal[index_bimodal] = pred_table_bimodal[index_bimodal] - 1;
				}

			}


			if (takenornot == 't')
				reality = 1;
			else if (takenornot == 'n')
				reality = 0;

			if (proceedtoBP == 0)
			{
				if (btb_miss == 1 && reality == 1)
					num_btb_missandtaken++;
			}

			if (proceedtoBP == 1 || (BTB_ASSOC == 0 && BTB_SIZE == 0))
			{
				if (prediction != reality)
					num_mispred_bp++;
			}
			proceedtoBP = 0;

			if (prediction != reality)
				misprediction++;

		}

		else if (!strcmp(argv[1], "gshare"))
		{
			isitgshare = 1;
			unsigned long initial = ((1 << M1) - 1) & (branch >> 2);
			unsigned long initial2 = (initial >> (M1 - N)) ^ (bhr);
			unsigned long sq = pow(2, N);
			unsigned long index = initial & ~((sq - 1) << (M1 - N));
			index = index | (initial2 << (M1 - N));
			if (BTB_SIZE != 0 && BTB_ASSOC != 0)
			{
				unsigned long tag = branch >> (bo + indexsize); //extract the tag from input address
				unsigned long btb_index = ((1 << indexsize) - 1) & (branch >> (bo));
				int r_hit_block = 0;
				int BTB_hit = 0;



				for (unsigned long j = 0; j < BTB_ASSOC; j++)
				{ //find if hit or miss and if hit find the hit block     
					if (btb_cache[btb_index][j].t_tag == tag)
					{
						BTB_hit = 1;
						r_hit_block = j;
					}
				}


				if (BTB_hit == 1)
				{ //hit case
					unsigned int oldLRU;
					oldLRU = btb_cache[btb_index][r_hit_block].t_LRU_counter;
					for (unsigned long i = 0; i < BTB_ASSOC; i++)
					{
						if (btb_cache[btb_index][i].t_LRU_counter < oldLRU) //update lru counters if the value of counters
							btb_cache[btb_index][i].t_LRU_counter += 1; //is less than hit block counter value
					}
					btb_cache[btb_index][r_hit_block].t_LRU_counter = 0; //put lru counter of hit block to zero
					proceedtoBP = 1;
				}




				if (BTB_hit == 0)
				{ //miss case 
					read_misses++;
					int invalid = 0;
					prediction = 0;
					btb_miss = 1;

					for (unsigned int vl = 0; vl < BTB_ASSOC; vl++)
					{
						if (btb_cache[btb_index][vl].t_valid_bit == 0)
						{ //check for number of invalids
							invalid++;
						}
					}


					if (invalid == 1)
					{ //case when only one block is empty
						for (unsigned int i = 0; i < BTB_ASSOC; i++)
						{
							if (btb_cache[btb_index][i].t_valid_bit == 0)
							{ //update tag,LRU_counter and valid_bit 
								for (unsigned int j = 0; j < BTB_ASSOC; j++)
								{
									btb_cache[btb_index][j].t_LRU_counter++;
								}

								btb_cache[btb_index][i].t_tag = tag;
								btb_cache[btb_index][i].t_LRU_counter = 0;
								btb_cache[btb_index][i].t_valid_bit = 1;
							}
						}
					}


					if (invalid == 0)
					{ //when all are filled, check highest LRU_count and update tag to that block
						unsigned int temp, highest_LRU;
						temp = btb_cache[btb_index][0].t_LRU_counter;
						highest_LRU = 0;
						for (unsigned int j = 0; j < BTB_ASSOC; j++)
						{
							if (btb_cache[btb_index][j].t_LRU_counter > temp)
							{
								temp = btb_cache[btb_index][j].t_LRU_counter;
								highest_LRU = j;
							}
						}

						for (unsigned int k = 0; k < BTB_ASSOC; k++)
						{
							btb_cache[btb_index][k].t_LRU_counter++;
						}

						btb_cache[btb_index][highest_LRU].t_tag = tag;
						btb_cache[btb_index][highest_LRU].t_LRU_counter = 0;
						btb_cache[btb_index][highest_LRU].t_dirty_bit = 0;

					}


					if (invalid != 0 && invalid != 1)
					{ //case when multiple blocks are empty     
						unsigned int currentHighest = 0, highest_LRU;
						for (unsigned int k = 0; k < BTB_ASSOC; k++)
						{
							if (btb_cache[btb_index][k].t_valid_bit == 0 && btb_cache[btb_index][k].t_LRU_counter > currentHighest)
							{ //check which are empty and find the one with highest lru count
								currentHighest = btb_cache[btb_index][k].t_LRU_counter;
								highest_LRU = k;
							}

						}

						//Update the LRU Counters
						for (unsigned int k = 0; k < BTB_ASSOC; k++)
						{
							if (btb_cache[btb_index][k].t_LRU_counter < btb_cache[btb_index][highest_LRU].t_LRU_counter)
								btb_cache[btb_index][k].t_LRU_counter++;
						}

						btb_cache[btb_index][highest_LRU].t_tag = tag; //update tag,lru counter and valid bit 
						btb_cache[btb_index][highest_LRU].t_LRU_counter = 0;
						btb_cache[btb_index][highest_LRU].t_valid_bit = 1;
					}

				}
			}

			if (proceedtoBP == 1 || (BTB_ASSOC == 0 && BTB_SIZE == 0))
			{
				num_pred_bp++;
				if (pred_table_gshare[index] >= 2)
					prediction = 1;
				else if (pred_table_gshare[index] < 2)
					prediction = 0;

				if (takenornot == 't')
				{
					if (pred_table_gshare[index] < 3)
						pred_table_gshare[index] = pred_table_gshare[index] + 1;
				}
				else if (takenornot == 'n')
				{
					if (pred_table_gshare[index] > 0)
						pred_table_gshare[index] = pred_table_gshare[index] - 1;
				}
			}

			if (takenornot == 't')
				reality = 1;
			else if (takenornot == 'n')
				reality = 0;

			if (proceedtoBP == 0)
			{
				if (btb_miss == 1 && reality == 1)
					num_btb_missandtaken++;
			}

			if (proceedtoBP == 1 || (BTB_ASSOC == 0 && BTB_SIZE == 0))
			{
				if (prediction != reality)
					num_mispred_bp++;
			}


			if (prediction != reality)
				misprediction++;

			if (proceedtoBP == 1 || (BTB_ASSOC == 0 && BTB_SIZE == 0))
			{
				bhr = bhr >> 1;
				unsigned long temp = N - 1;
				bhr = bhr & ~(1 << temp);
				bhr = bhr | (reality << temp);
				proceedtoBP = 0;
			}
		}

		else if (!strcmp(argv[1], "hybrid"))
		{
			isithybrid = 1;
			unsigned long index = ((1 << K) - 1) & (branch >> 2);

			unsigned long index_bimodal = ((1 << M2) - 1) & (branch >> 2);

			if (pred_table_bimodal[index_bimodal] >= 2)
				prediction_bimodal = 1;
			else if (pred_table_bimodal[index_bimodal] < 2)
				prediction_bimodal = 0;

			unsigned long initial_gshare = ((1 << M1) - 1) & (branch >> 2);
			unsigned long initial2_gshare = (initial_gshare >> (M1 - N)) ^ (bhr);
			unsigned long sq = pow(2, N);
			unsigned long index_gshare = initial_gshare & ~((sq - 1) << (M1 - N));
			index_gshare = index_gshare | (initial2_gshare << (M1 - N));

			if (pred_table_gshare[index_gshare] >= 2)
				prediction_gshare = 1;
			else if (pred_table_gshare[index_gshare] < 2)
				prediction_gshare = 0;

			if (chooser_table[index] >= 2)
				prediction = prediction_gshare;
			else if (chooser_table[index] < 2)
				prediction = prediction_bimodal;


			if (takenornot == 't')
			{
				reality = 1;
				if (chooser_table[index] < 2)
				{
					if (pred_table_bimodal[index_bimodal] != 3)
						pred_table_bimodal[index_bimodal] = pred_table_bimodal[index_bimodal] + 1;
				}

				else if (chooser_table[index] >= 2)
				{
					if (pred_table_gshare[index_gshare] < 3)
						pred_table_gshare[index_gshare] = pred_table_gshare[index_gshare] + 1;
				}
			}

			else if (takenornot == 'n')
			{
				reality = 0;
				if (chooser_table[index] < 2)
				{
					if (pred_table_bimodal[index_bimodal] != 0)
						pred_table_bimodal[index_bimodal] = pred_table_bimodal[index_bimodal] - 1;
				}
				else if (chooser_table[index] >= 2)
				{
					if (pred_table_gshare[index_gshare] != 0)
						pred_table_gshare[index_gshare] = pred_table_gshare[index_gshare] - 1;
				}
			}
			if (prediction != reality)
				misprediction++;

			bhr = bhr >> 1;
			unsigned long temp = N - 1;
			bhr = bhr & ~(1 << temp);
			bhr = bhr | (reality << temp);

			if (reality == prediction_bimodal && reality != prediction_gshare && chooser_table[index] > 0)
				chooser_table[index] = chooser_table[index] - 1;
			else if (reality != prediction_bimodal && reality == prediction_gshare && chooser_table[index] < 3)
				chooser_table[index] = chooser_table[index] + 1;
		}



	}
	double mispred_rate = 0.0 f;
	mispred_rate = (misprediction * 100.0 / num_predictions);

	if (BTB_SIZE == 0 && BTB_ASSOC == 0)
	{
		if (isitbimodal == 1)
		{
			cout << "COMMAND" << endl;
			cout << "./sim " << predictor << " " << M2 << " " << BTB_SIZE << " " << BTB_ASSOC << " " << trace_file << endl;
			cout << "OUTPUT" << endl;
			cout << "number of predictions: " << num_predictions << endl;
			std::cout << std::fixed;
			std::cout << std::setprecision(2);
			cout << "number of mispredictions: " << misprediction << endl;
			cout << "misprediction rate: " << mispred_rate << "%" << endl;
			cout << "FINAL BIMODAL CONTENTS" << endl;
			for (unsigned long m = 0; m < pow(2, M2); m++)
				cout << m << " " << pred_table_bimodal[m] << endl;
		}
		else if (isitgshare == 1)
		{
			cout << "COMMAND" << endl;
			cout << "./sim " << predictor << "" << M1 << " " << N << " " << BTB_SIZE << " " << BTB_ASSOC << " " << trace_file << endl;
			cout << "OUTPUT" << endl;
			cout << "number of predictions: " << num_predictions << endl;
			std::cout << std::fixed;
			std::cout << std::setprecision(2);
			cout << "number of mispredictions: " << misprediction << endl;
			cout << "misprediction rate: " << mispred_rate << "%" << endl;
			cout << "FINAL GSHARE CONTENTS" << endl;
			for (unsigned long m = 0; m < pow(2, M1); m++)
				cout << m << " " << pred_table_gshare[m] << endl;
		}
		else if (isithybrid == 1)
		{
			cout << "COMMAND" << endl;
			cout << "./sim " << predictor << " " << K << " " << M1 << " " << N << " " << M2 << " " << BTB_SIZE << " " << BTB_ASSOC << " " << trace_file << endl;
			cout << "OUTPUT" << endl;
			cout << "number of predictions: " << num_predictions << endl;
			std::cout << std::fixed;
			std::cout << std::setprecision(2);
			cout << "number of mispredictions: " << misprediction << endl;
			cout << "misprediction rate: " << mispred_rate << "%" << endl;
			cout << "FINAL CHOOSER CONTENTS" << endl;
			for (unsigned long m = 0; m < pow(2, K); m++)
				cout << m << " " << chooser_table[m] << endl;
			cout << "FINAL GSHARE CONTENTS" << endl;
			for (unsigned long m = 0; m < pow(2, M1); m++)
				cout << m << " " << pred_table_gshare[m] << endl;
			cout << "FINAL BIMODAL CONTENTS" << endl;
			for (unsigned long m = 0; m < pow(2, M2); m++)
				cout << m << " " << pred_table_bimodal[m] << endl;
		}
	}
	else if (BTB_SIZE != 0 && BTB_ASSOC != 0)
	{
		if (isitbimodal == 1)
		{
			cout << "COMMAND" << endl;
			cout << "./sim " << predictor << " " << M2 << " " << BTB_SIZE << " " << BTB_ASSOC << " " << trace_file << endl;
			cout << "OUTPUT" << endl;
			cout << "size of BTB: " << BTB_SIZE << endl;
			cout << "number of branches: " << num_predictions << endl;
			cout << "number of predictions from branch predictor: " << num_pred_bp << endl;
			cout << "number of mispredictions from branch predictor: " << num_mispred_bp << endl;
			cout << "number of branches miss in BTB and taken: " << num_btb_missandtaken << endl;
			std::cout << std::fixed;
			std::cout << std::setprecision(2);
			cout << "total mispredictions: " << misprediction << endl;
			cout << "misprediction rate: " << mispred_rate << "%" << endl;
			cout << "FINAL BTB CONTENTS" << endl;
			for (unsigned int i = 0; i < sets; i++)
			{
				cout << "set  " << dec << i << ": ";
				for (unsigned int j = 0; j < BTB_ASSOC; j++)
				{
					for (unsigned int k = 0; k < BTB_ASSOC; k++)
					{
						if (btb_cache[i][k].t_LRU_counter == j)
						{
							cout << " " << hex << btb_cache[i][k].t_tag << " ";
						}
					}
				}
				cout << endl;
			}
			cout << endl;
			cout << "FINAL BIMODAL CONTENTS" << endl;
			for (unsigned long m = 0; m < pow(2, M2); m++)
				cout << dec << m << " " << pred_table_bimodal[m] << endl;

		}
		else if (isitgshare == 1)
		{
			cout << "COMMAND" << endl;
			cout << "./sim " << predictor << " " << M1 << " " << N << " " << BTB_SIZE << " " << BTB_ASSOC << " " << trace_file << endl;
			cout << "OUTPUT" << endl;
			cout << "size of BTB: " << BTB_SIZE << endl;
			cout << "number of branches: " << num_predictions << endl;
			cout << "number of predictions from branch predictor: " << num_pred_bp << endl;
			cout << "number of mispredictions from branch predictor: " << num_mispred_bp << endl;
			cout << "number of branches miss in BTB and taken: " << num_btb_missandtaken << endl;
			std::cout << std::fixed;
			std::cout << std::setprecision(2);
			cout << "total mispredictions: " << misprediction << endl;
			cout << "misprediction rate: " << mispred_rate << "%" << endl;
			cout << "FINAL BTB CONTENTS" << endl;
			for (unsigned int i = 0; i < sets; i++)
			{
				cout << "set  " << dec << i << ": ";
				for (unsigned int j = 0; j < BTB_ASSOC; j++)
				{
					for (unsigned int k = 0; k < BTB_ASSOC; k++)
					{
						if (btb_cache[i][k].t_LRU_counter == j)
						{
							if (btb_cache[i][k].t_tag != 0)
								cout << " " << hex << btb_cache[i][k].t_tag << " ";
							else if (btb_cache[i][k].t_tag == 0)
								cout << " ";
						}
					}
				}
				cout << endl;
			}
			cout << endl;
			cout << "FINAL GSHARE CONTENTS" << endl;
			for (unsigned long m = 0; m < pow(2, M1); m++)
				cout << dec << m << " " << pred_table_gshare[m] << endl;
		}
	}
	return 0;
}
