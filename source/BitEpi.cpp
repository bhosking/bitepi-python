
#ifdef _MSC_VER
	#define _CRT_SECURE_NO_WARNINGS
	typedef int pthread_mutex_t;
	int pthread_mutex_trylock(pthread_mutex_t *x)
	{
		if (*x == 0)
		{
			*x = 1;
			return 0;
		}
		else
			return 1;
	}
	void pthread_mutex_init(pthread_mutex_t *x, int *y)
	{
		*x = 0;
	}
	typedef int pthread_t;
	typedef int pthread_attr_t;
	void pthread_create(pthread_t *thread, const pthread_attr_t *attr, void *(*start_routine) (void *), void *arg)
	{
		(*start_routine)(arg);
	}
	void pthread_join(pthread_t thread, void **retval)
	{
		return;
	}
#else
	#include "pthread.h"
#endif

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "time.h"
#include "math.h"
#include "csvparser.h"

//#define PTEST

#ifdef PTEST
	clock_t elapse[100];
#endif


typedef unsigned char uint8;
typedef unsigned short int uint16;
typedef unsigned int uint32;
typedef unsigned long long int uint64;
typedef int int32;

typedef unsigned int varIdx;
typedef unsigned short int sampleIdx; // this type used in contingency table. this table should be kept in cache. so choose the smallest possible type here
typedef unsigned long long int word; // for parallel processing

const uint8 cti[81] = {0,1,2,4,5,6,8,9,10,16,17,18,20,21,22,24,25,26,32,33,34,36,37,38,40,41,42,64,65,66,68,69,70,72,73,74,80,81,82,84,85,86,88,89,90,96,97,98,100,101,102,104,105,106,128,129,130,132,133,134,136,137,138,144,145,146,148,149,150,152,153,154,160,161,162,164,165,166,168,169,170 };
const uint32 byte_in_word = sizeof(word);

#define MAX_ORDER 4

#define P2(X) (X*X)
#define P3(X) (X*X*X)
#define P4(X) (X*X*X*X)

#define ERROR(X) {printf("*** ERROR: %s (line:%u - File %s)\n", X, __LINE__, __FILE__); exit(0);}
#define NULL_CHECK(X) {if(!X) {printf("*** ERROR: %s is null (line:%u - File %s)\n", #X, __LINE__, __FILE__); exit(0);}}

double ***tripletPurity;
double **PairPurity;
double *SnpPurity;

union WordByte
{
	word w;
	uint8 b[8];
};

struct ARGS
{
	bool computeP[MAX_ORDER]; // [N] should we compute purity of order of N
	bool printP[MAX_ORDER];   // [N] should we report purity of order of N if it meets the threshold p[N]
	bool saveP[MAX_ORDER];    // [N] should we save purity of order of N to compute IG of order N+1
	bool computeIG[MAX_ORDER];// [N] should we compute IG of order of N
	bool printIG[MAX_ORDER];  // [N] should we report IG of order of N in it meets the threshold ig[N]
	bool bestIG;			  // [N] should we compute the best IG

	double p[MAX_ORDER];
	double ig[MAX_ORDER];

	char input[1024];
	char output[1024];
	uint32 numThreads;
	uint32 order;

	bool sort;

	ARGS()
	{
		memset(this, 0, sizeof(ARGS));
		numThreads = 1;
		order = 1;
	}

	~ARGS()
	{
	}

	void PrintHelp(char* exec)
	{
		printf(" -i		Input CSV file\n");
		printf("		* First row includes labels 1 and 0 for case and controls\n");
		printf("		* First column includes SNP ids\n");
		printf("		* First entry (first col and first row) is ignored\n");
		printf("		* All other entry can be 0, 1 or 2 (HomRef, Het and HomVar genotype respectively)\n");
		
		printf(" -o		Output prefix\n");

		printf(" -sort		Sort output files by Purity and Information-Gained\n");
		
		printf(" -t		number of threads\n");

		printf(" -bestIG	find the best interactions for each SNP (will disregards below options)\n");

		printf(" -p1 [thr]	Compute purity for 1-SNP (SNP).\n");
		printf(" -p2 [thr]	Compute purity for 2-SNP (Pair).\n");
		printf(" -p3 [thr]	Compute purity for 3-SNP (Triplet).\n");
		printf(" -p4 [thr]	Compute purity for 4-SNP (Quadlet).\n");

		printf(" -ig1 [thr]	Compute Information-Gained (IG) for 1-SNP (SNP).\n");
		printf(" -ig2 [thr]	Compute Information-Gained (IG) for 2-SNP (Pair).\n");
		printf(" -ig3 [thr]	Compute Information-Gained (IG) for 3-SNP (Triplet).\n");
		printf(" -ig4 [thr]	Compute Information-Gained (IG) for 4-SNP (Quadlet).\n");

		printf("* thr is threshold and is optional. If you dont pass thr it computes the metric but it does not report anything (performance testing).\n");
		printf("* 0<thr<1.\n");
		printf("* if you want all interactions set thr to 0.\n");
		printf("* if you set ig(n) it will compute p(n-1) as it is needed.\n");

		printf("\n====================================\n");
		
		ERROR("Please enter valid arguments");
		return;
	}

	void Parse(int argc, char* argv[])
	{
		double d = -1;
		uint32 o = 0;
		char str[100];
		bool next = false;

		for (uint32 i = 1; i < argc; i++)
		{
			next = false;

			// check if purity flag is passed for any order
			for (uint32 o = 0; o < MAX_ORDER; o++)
			{
				sprintf(str, "-p%u", o+1);
				if (!strcmp(argv[i], str))
				{
					// set the computep flag
					computeP[o] = true;

					// check if there is any threashold argument to this option
					if ((i + 1) != argc)
					{
						if (argv[i + 1][0] != '-')
						{
							d = atof(argv[i + 1]);
							if ((d == 0) && (argv[i + 1][0] != '0'))
								PrintHelp(argv[0]);
							// set the threshold and print flag for purity
							p[o] = d;
							printP[o] = true;
							i++;
						}
					}
					next = true;
					break;
				}
			}
			// if this option is processed move to the next option
			if (next) continue;

			// check if ig flag is passed for any order
			for (uint32 o = 0; o < MAX_ORDER; o++)
			{
				sprintf(str, "-ig%u", o + 1);
				if (!strcmp(argv[i], str))
				{
					// set the computep flag and purity flag of previous order
					computeP[o] = computeIG[o] = true;
					if(o>0)
						computeP[o-1] = saveP[o-1] = true;

					// check if there is any threashold argument to this option
					if ((i + 1) != argc)
					{
						if (argv[i + 1][0] != '-')
						{
							d = atof(argv[i + 1]);
							if ((d == 0) && (argv[i + 1][0] != '0'))
								PrintHelp(argv[0]);
							// set the threshold and print flag for ig
							ig[o] = d;
							printIG[o] = true;
							i++;
						}
					}
					next = true;
					break;
				}
			}
			if (next) continue;

			// read input file name
			if (!strcmp(argv[i], "-i"))
			{
				if ((i + 1) == argc)
					PrintHelp(argv[0]);

				if (argv[i + 1][0] != '-')
					strcpy(input, argv[i+1]);
				else
					PrintHelp(argv[0]);
				i++;
				continue;
			}

			// read output file prefix
			if (!strcmp(argv[i], "-o"))
			{
				if ((i + 1) == argc)
					PrintHelp(argv[0]);

				if (argv[i + 1][0] != '-')
					strcpy(output, argv[i + 1]);
				else
					PrintHelp(argv[0]);
				i++;
				continue;
			}

			// read number of threads
			if (!strcmp(argv[i], "-t"))
			{
				if ((i + 1) == argc)
					PrintHelp(argv[0]);

				if (argv[i + 1][0] != '-')
				{
					numThreads = atoi(argv[i + 1]);
					if ((numThreads == 0) && (argv[i + 1][0] != '0'))
						PrintHelp(argv[0]);
				}
				else
					PrintHelp(argv[0]);
				i++;
				continue;
			}

			// read bestIG flag
			if (!strcmp(argv[i], "-bestIG"))
			{
				bestIG = true;
				continue;
			}

			// read bestIG flag
			if (!strcmp(argv[i], "-sort"))
			{
				sort = true;
				continue;
			}

			printf("\n***ERR*** invalid option %s\n", argv[i]);
			PrintHelp(argv[0]);
		}

		// check arguments
		if (strlen(input) == 0 || strlen(output) == 0)
			PrintHelp(argv[0]);

		// apply bestIG
		if(bestIG)
		for (uint32 o = 0; o < MAX_ORDER; o++)
		{
			computeP[o] = saveP[o] = computeIG[o] = true;
			printP[o] = printIG[o] = false;
		}

		if (computeP[0]) order = 1;
		if (computeP[1]) order = 2;
		if (computeP[2]) order = 3;
		if (computeP[3]) order = 4;
	}

	void Print()
	{
		printf("\n -i		%s", input);
		printf("\n -o		%s", output);
		printf("\n -t		%u", numThreads);
		printf("\n order		%u", order);
		printf("\n bestIG		%s", bestIG ? "true" : "false");
		printf("\n sort			%s", sort ? "true" : "false");
	
		for (uint32 o = 0; o < MAX_ORDER; o++)
		{
			printf("\n -p%u		%f", o, p[o]);
			printf("\n -ig%u		%f", o, ig[o]);
			printf("\n computeP[%u]	%s", o, computeP[o] ? "true" : "false");
			printf("\n printP[%u]	%s", o, printP[o] ? "true" : "false");
			printf("\n saveP[%u]	%s", o, saveP[o] ? "true" : "false");
			printf("\n computeIG[%u]	%s", o, computeIG[o] ? "true" : "false");
			printf("\n printIG[%u]	%s", o, printIG[o]?"true":"false");
		}
	}
};

void AllocatePurity(varIdx n, ARGS args)
{
	if (args.saveP[0])
	{
		SnpPurity = new double[n];
		NULL_CHECK(SnpPurity);
	}

	if (args.saveP[1])
	{
		PairPurity = new double *[n];
		NULL_CHECK(PairPurity);
		for (varIdx i = 0; i < n; i++)
		{
			PairPurity[i] = new double[n];
			NULL_CHECK(PairPurity[i]);
		}
	}

	if (args.saveP[2])
	{
		tripletPurity = new double **[n];
		NULL_CHECK(tripletPurity);
		
		for (varIdx i = 0; i < n; i++)
		{
			tripletPurity[i] = new double *[n];
			NULL_CHECK(tripletPurity[i]);

			for (varIdx j = 0; j < n; j++)
			{
				tripletPurity[i][j] = new double[n];
				NULL_CHECK(tripletPurity[i][j]);
			}

		}
	}
}

void FreePurity(varIdx n, ARGS args)
{
	if (args.saveP[0])
	{
		delete[] SnpPurity;
	}

	if (args.saveP[1])
	{
		for (varIdx i = 0; i < n; i++)
		{
			delete[] PairPurity[i];
		}
		delete[] PairPurity;
	}

	if (args.saveP[2])
	{
		for (varIdx i = 0; i < n; i++)
		{
			for (varIdx j = 0; j < n; j++)
			{
				delete[] tripletPurity[i][j];
			}
			delete[] tripletPurity[i];
		}
		delete[] tripletPurity;
	}
}

struct InformationGained
{
public:
	double purity[4];
	double ig[4];
	varIdx pair[2];
	varIdx triplet[3];
	varIdx quadlet[4];

	void toCSV(FILE *csv, varIdx id, char **names)
	{
		fprintf(csv, "%s,", names[id]);
		fprintf(csv, "%f,%f,%f,%f,", purity[0], purity[1], purity[2], purity[3]);
		fprintf(csv, "%f,%f,%f,%f,", ig[0], ig[1], ig[2], ig[3]);
		
		if (id == pair[0])
			fprintf(csv, "%s,", names[pair[1]]);
		else
			fprintf(csv, "%s,", names[pair[0]]);

		if (id == triplet[0])
			fprintf(csv, "%s,%s,", names[triplet[1]], names[triplet[2]]);
		else if (id == triplet[1])
			fprintf(csv, "%s,%s,", names[triplet[0]], names[triplet[2]]);
		else
			fprintf(csv, "%s,%s,", names[triplet[0]], names[triplet[1]]);

		if (id == quadlet[0])
			fprintf(csv, "%s,%s,%s\n", names[quadlet[1]], names[quadlet[2]], names[quadlet[3]]);
		else if (id == quadlet[1])
			fprintf(csv, "%s,%s,%s\n", names[quadlet[0]], names[quadlet[2]], names[quadlet[3]]);
		else if (id == quadlet[2])
			fprintf(csv, "%s,%s,%s\n", names[quadlet[0]], names[quadlet[1]], names[quadlet[3]]);
		else
			fprintf(csv, "%s,%s,%s\n", names[quadlet[0]], names[quadlet[1]], names[quadlet[2]]);
	}

	void Max(const InformationGained &o)
	{
		if (o.ig[0] > ig[0])
		{
			ig[0] = o.ig[0];
			purity[0] = o.purity[0];
		}

		if (o.ig[1] > ig[1])
		{
			ig[1] = o.ig[1];
			purity[1] = o.purity[1];
			pair[0] = o.pair[0];
			pair[1] = o.pair[1];
		}

		if (o.ig[2] > ig[2])
		{
			ig[2] = o.ig[2];
			purity[2] = o.purity[2];
			triplet[0] = o.triplet[0];
			triplet[1] = o.triplet[1];
			triplet[2] = o.triplet[2];
		}

		if (o.ig[3] > ig[3])
		{
			ig[3] = o.ig[3];
			purity[3] = o.purity[3];
			quadlet[0] = o.quadlet[0];
			quadlet[1] = o.quadlet[1];
			quadlet[2] = o.quadlet[2];
			quadlet[3] = o.quadlet[3];
		}
	}
};

class Result
{
public:
	varIdx numVariable;
	InformationGained *res;

	~Result()
	{
		delete[] res;
	}

	void toCSV(char *fn, char **names)
	{
		FILE *csv = fopen(fn, "w");
		NULL_CHECK(csv);

		fprintf(csv, "SNP,SNP_P,PAIR_P,TRIPLET_P,QUADLET_P,SNP_IG,PAIR_IG,TRIPLET_IG,QUADLET_IG,PAIR,TRIPLET_1,TRIPLET_2,QUADLET_1,QUADLET_2,QUADLET_3\n");

		for (varIdx i = 0; i < numVariable; i++)
			res[i].toCSV(csv, i, names);

		fclose(csv);
	}

	void Init(varIdx nv)
	{
		numVariable = nv;
		res = new InformationGained[numVariable];
		NULL_CHECK(res);
		memset(res, 0, numVariable * sizeof(InformationGained));
	}

	void Max(const Result &o)
	{
		for (varIdx i = 0; i < numVariable; i++)
		{
			res[i].Max(o.res[i]);
		}
	}

	void Max_1(double ig, double purity, varIdx *idx)
	{
		if (ig > res[idx[0]].ig[0])
		{
			res[idx[0]].ig[0] = ig;
			res[idx[0]].purity[0] = purity;
		}
	}

	void Max_2(double ig, double purity, varIdx *idx)
	{
		InformationGained *r;
		for (uint32 i = 0; i < 2; i++)
		{
			r = &res[idx[i]];
			if (ig > r->ig[1])
			{
				r->ig[1] = ig;
				r->purity[1] = purity;
				r->pair[0] = idx[0];
				r->pair[1] = idx[1];
			}
		}
	}

	void Max_3(double ig, double purity, varIdx *idx)
	{
		InformationGained *r;
		for (uint32 i = 0; i < 3; i++)
		{
			r = &res[idx[i]];
			if (ig > r->ig[2])
			{
				r->ig[2] = ig;
				r->purity[2] = purity;
				r->triplet[0] = idx[0];
				r->triplet[1] = idx[1];
				r->triplet[2] = idx[2];
			}
		}
	}

	void Max_4(double ig, double purity, varIdx *idx)
	{
		InformationGained *r;
		for (uint32 i = 0; i < 4; i++)
		{
			r = &res[idx[i]];
			if (ig > r->ig[3])
			{
				r->ig[3] = ig;
				r->purity[3] = purity;
				r->quadlet[0] = idx[0];
				r->quadlet[1] = idx[1];
				r->quadlet[2] = idx[2];
				r->quadlet[3] = idx[3];
			}
		}
	}

};

class Dataset
{
	// contigency table index translation
	// note that 3 or 0b11 is not a valid genotype and should not be considered in Gini Computation.
	uint64 CaseIndex(varIdx v, sampleIdx s) { return ((v * numByteCase) + s); }	// get the byte index of sample in data
	uint64 CtrlIndex(varIdx v, sampleIdx s) { return ((v * numByteCtrl) + s); }	// get the byte index of sample in data

public:

	uint32 order;

	uint32 *labels;

	sampleIdx numSample;

	sampleIdx numCase;
	sampleIdx numCtrl;

	uint32 numWordCase;// number of machine word used to store Case data (each sample is a byte)
	uint32 numWordCtrl;// number of machine word used to store Ctrl data (each sample is a byte)

	uint32 numByteCase; // numWordCase * sizeof(word)
	uint32 numByteCtrl; // numWordCtrl * sizeof(word)

	uint8 *byteCase[MAX_ORDER]; // byte pointer to store genotype data and their shifted version
	uint8 *byteCtrl[MAX_ORDER]; // byte pointer to store genotype data and their shifted version

	word *wordCase[MAX_ORDER]; // the wrod pointer to byteCase
	word *wordCtrl[MAX_ORDER]; // the wrod pointer to byteCtrl

	uint32 numLine;
	varIdx numVariable;
	char **nameVariable;

	double setPurity; // purity of the original set

	sampleIdx *contingency_table; // should be small enough to remain in cache

	Result *results;

	void FreeMemory(ARGS args)
	{
		delete[] labels;

		for (uint32 i = 0; i < numVariable; i++)
			delete[] nameVariable[i];
		delete nameVariable;

		for (uint32 i = 0; i < order; i++)
		{
			delete[] wordCase[i];
			delete[] wordCtrl[i];
		}
		
		if(args.bestIG)
			delete[] results;
	}

	// Count number of line in a file to see how many variable exits.
	uint32 LineCount(const char *fn)
	{
		printf("\nCounting lines in %s.", fn);

		FILE *f;
		uint32 lines = 0;

		f = fopen(fn, "r");
		NULL_CHECK(f)

			char ch;
		for (ch = getc(f); ch != EOF; ch = getc(f)) if (ch == '\n') lines = lines + 1;

		fclose(f);

		printf("\n%s has %u lines\n", fn, lines);
		return lines;
	}

	// This function read data from file
	void ReadDataset(const char *fn)
	{
		printf("\nloading dataset %s", fn);

		numLine = LineCount(fn);
		nameVariable = new char*[numLine - 1];
		NULL_CHECK(nameVariable);

		CsvParser *csvparser = CsvParser_new(fn, ",", 1);
		CsvRow *row;

		const CsvRow *header = CsvParser_getHeader(csvparser);
		NULL_CHECK(header);

		const char **headerFields = CsvParser_getFields(header);

		numSample = CsvParser_getNumFields(header) - 1;
		if (numSample >= pow(2, sizeof(sampleIdx) * 8))
			ERROR("Change sampleIdx type to support number of samples exist in dataset");

		labels = new uint32[numSample];
		NULL_CHECK(labels);

		numCase = 0;
		numCtrl = 0;
		for (sampleIdx i = 0; i < numSample; i++)
		{
			sscanf(headerFields[i + 1], "%u", &labels[i]);
			if (labels[i] > 1) ERROR("Class shold be 0 or 1");

			if (labels[i])
				numCase++;
			else
				numCtrl++;
		}

		// find number of word and byte per variable in Case and Ctrl
		numWordCase = numCase / byte_in_word;
		numWordCtrl = numCtrl / byte_in_word;

		if (numCase % byte_in_word) numWordCase++;
		if (numCtrl % byte_in_word) numWordCtrl++;

		numByteCase = numWordCase * byte_in_word;
		numByteCtrl = numWordCtrl * byte_in_word;

		// allocate memory
		wordCase[0] = new word[numLine * numWordCase];
		wordCtrl[0] = new word[numLine * numWordCtrl];

		NULL_CHECK(wordCase[0]);
		NULL_CHECK(wordCtrl[0]);

		// convert to byte address
		byteCase[0] = (uint8*)wordCase[0];
		byteCtrl[0] = (uint8*)wordCtrl[0];

		numVariable = 0;
		while ((row = CsvParser_getRow(csvparser)))
		{
			const char **rowFields = CsvParser_getFields(row);

			if (CsvParser_getNumFields(row) != (numSample + 1))
				ERROR("Number of fields does not match the first line in the file");

			nameVariable[numVariable] = new char[strlen(rowFields[0]) + 1];
			strcpy(nameVariable[numVariable], rowFields[0]);

			uint32 idxCase = 0;
			uint32 idxCtrl = 0;

			for (sampleIdx i = 0; i < numSample; i++)
			{
				uint32 gt;
				sscanf(rowFields[i + 1], "%u", &gt);
				if (gt > 2) ERROR("Values shold be 0 or 1 or 2");
				if (labels[i])
				{
					byteCase[0][CaseIndex(numVariable, idxCase)] = (uint8)gt;
					idxCase++;
				}
				else
				{
					byteCtrl[0][CtrlIndex(numVariable, idxCtrl)] = (uint8)gt;
					idxCtrl++;
				}
			}
			CsvParser_destroy_row(row);

			numVariable++;
		}

		CsvParser_destroy(csvparser);

		return;
	}

	// This function write data from file (to test ReadDataset function)
	void WriteDataset(const char *fn)
	{
		printf("\nWrite dataset to %s", fn);

		FILE *f = fopen(fn, "w");
		NULL_CHECK(f);

		fprintf(f, "Var\\Class");

		for (uint32 i = 0; i < numSample; i++)
		{
			fprintf(f, ",%u", labels[i]);
		}
		fprintf(f, "\n");

		for (uint32 j = 0; j < numVariable; j++)
		{
			fprintf(f, "%s", nameVariable[j]);

			uint32 idxCase = 0;
			uint32 idxCtrl = 0;

			for (uint32 i = 0; i < numSample; i++)
			{
				if (labels[i])
				{
					fprintf(f, ",%u", byteCase[0][CaseIndex(j, idxCase)]);
					idxCase++;
				}
				else
				{
					fprintf(f, ",%u", byteCtrl[0][CtrlIndex(j, idxCtrl)]);
					idxCtrl++;
				}
			}
			fprintf(f, "\n");
		}
		fclose(f);
		return;
	}

	void ComputeSetPurity()
	{
		numCase = 0;
		numCtrl = 0;
		for (sampleIdx i = 0; i < numSample; i++)
			if (labels[i]) numCase++; else numCtrl++;
		setPurity = P2((double)numCase / numSample) + P2((double)numCtrl / numSample);
		printf("\nSet Purity of the dataset is %f", setPurity);
	}

	void Shift()
	{
		for (uint32 d = 1; d < order; d++)
		{
			// allocate memory
			wordCase[d] = new word[numVariable * numWordCase];
			wordCtrl[d] = new word[numVariable * numWordCtrl];
			NULL_CHECK(wordCase[d]);
			NULL_CHECK(wordCtrl[d]);
			// convert to byte address
			byteCase[d] = (uint8*)wordCase[d];
			byteCtrl[d] = (uint8*)wordCtrl[d];

			// shoft and copy
			for (uint32 i = 0; i < (numVariable * numWordCase); i++)
				wordCase[d][i] = wordCase[d - 1][i] << 2;
			for (uint32 i = 0; i < (numVariable * numWordCtrl); i++)
				wordCtrl[d][i] = wordCtrl[d - 1][i] << 2;

			printf("\nShift dataset by %u bits compeleted", d * 2);
		}
	}

	void Init(ARGS args)
	{
		order = args.order;
		if (args.bestIG)
		{
			results = new Result[args.numThreads]; // number of threads
			for (uint32 i = 0; i < args.numThreads; i++)
				results[i].Init(numVariable);
		}
		ComputeSetPurity();
		Shift();
	}

	word *GetVarCase(uint32 o, varIdx vi)
	{
		return &wordCase[o][vi*numWordCase];
	}

	word *GetVarCtrl(uint32 o, varIdx vi)
	{
		return &wordCtrl[o][vi*numWordCtrl];
	}
};

struct ThreadData
{
	void *epiStat; // epi class
	uint32 id; // thread index
};

class EpiStat
{
public:
	Dataset *dataset;

	ARGS args;
	uint32 threadIdx;

	void *(*threadFunction[4]) (void *);

	// These 4 items must be allocated by each thread separately
	word *epiCaseWord[3];
	word *epiCtrlWord[3];

	sampleIdx *contingencyCase;
	sampleIdx *contingencyCtrl;

	FILE **topPfile;
	FILE **topIGfile;

	void OpenFiles(uint32 order)
	{
		topPfile = new FILE*[args.numThreads];
		NULL_CHECK(topPfile);

		topIGfile = new FILE*[args.numThreads];
		NULL_CHECK(topIGfile);

		char* fn = new char[strlen(args.output) + 20];
		NULL_CHECK(fn);
		for (uint32 t = 0; t < args.numThreads; t++)
		{
			if (args.printP[order])
			{
				sprintf(fn, "%s.Purity.%u.%u.csv", args.output, order, t);
				topPfile[t] = fopen(fn, "w");
				NULL_CHECK(topPfile)
			}
			if (args.printIG[order])
			{
				sprintf(fn, "%s.IG.%u.%u.csv", args.output, order, t);
				topIGfile[t] = fopen(fn, "w");
				NULL_CHECK(topIGfile)
			}
		}
		delete[]fn;
	}

	void CloseFiles(uint32 order)
	{
		for (uint32 t = 0; t < args.numThreads; t++)
		{
			if (args.printP[order])
			{
				fclose(topPfile[t]);
			}
			if (args.printIG[order])
			{
				fclose(topIGfile[t]);
			}
		}
	}

	EpiStat()
	{
		return;
	}

	~EpiStat()
	{
		
	}
	
	EpiStat(EpiStat *ref)
	{
		memcpy(this, ref, sizeof(EpiStat));
	}

	void Init(Dataset *d, ARGS a, void *(*tf1) (void *), void *(*tf2) (void *), void *(*tf3) (void *), void *(*tf4) (void *))
	{
		dataset = d;
		args = a;
		threadFunction[0] = tf1;
		threadFunction[1] = tf2;
		threadFunction[2] = tf3;
		threadFunction[3] = tf4;
	}

	void AllocateThreadMemory()
	{
		for (uint32 i = 0; i < MAX_ORDER-1; i++)
		{
			epiCaseWord[i] = new word[dataset->numWordCase];
			epiCtrlWord[i] = new word[dataset->numWordCtrl];

			NULL_CHECK(epiCaseWord[i]);
			NULL_CHECK(epiCtrlWord[i]);
		}

		contingencyCase = new sampleIdx[(uint32)pow(2, MAX_ORDER * 2)];
		contingencyCtrl = new sampleIdx[(uint32)pow(2, MAX_ORDER * 2)];

		NULL_CHECK(contingencyCase);
		NULL_CHECK(contingencyCtrl);
	}

	void FreeThreadMemory()
	{
		for (uint32 i = 0; i < MAX_ORDER-1; i++)
		{
			delete[] epiCaseWord[i];
			delete[] epiCtrlWord[i];
		}

		delete[] contingencyCase;
		delete[] contingencyCtrl;
	}

	void OR_1(varIdx idx)
	{
		const uint32 OIDX = 0; // SNP
		word *caseData = dataset->GetVarCase(OIDX, idx);
		word *ctrlData = dataset->GetVarCtrl(OIDX, idx);

		for (uint32 i = 0; i < dataset->numWordCase; i++)
		{
			epiCaseWord[OIDX][i] = caseData[i];
		}

		for (uint32 i = 0; i < dataset->numWordCtrl; i++)
		{
			epiCtrlWord[OIDX][i] = ctrlData[i];
		}
	}

	void OR_2(varIdx idx)
	{
		const uint32 OIDX = 1; // Pair
		word *caseData = dataset->GetVarCase(OIDX, idx);
		word *ctrlData = dataset->GetVarCtrl(OIDX, idx);

		for (uint32 i = 0; i < dataset->numWordCase; i++)
		{
			epiCaseWord[OIDX][i] = epiCaseWord[OIDX - 1][i] | caseData[i];
		}

		for (uint32 i = 0; i < dataset->numWordCtrl; i++)
		{
			epiCtrlWord[OIDX][i] = epiCtrlWord[OIDX - 1][i] | ctrlData[i];
		}
	}

	void OR_3(varIdx idx)
	{
		const uint32 OIDX = 2; // Triplet
		word *caseData = dataset->GetVarCase(OIDX, idx);
		word *ctrlData = dataset->GetVarCtrl(OIDX, idx);

		for (uint32 i = 0; i < dataset->numWordCase; i++)
		{
			epiCaseWord[OIDX][i] = epiCaseWord[OIDX - 1][i] | caseData[i];
		}

		for (uint32 i = 0; i < dataset->numWordCtrl; i++)
		{
			epiCtrlWord[OIDX][i] = epiCtrlWord[OIDX - 1][i] | ctrlData[i];
		}
	}

	void OR_1x(varIdx idx)
	{
		const uint32 OIDX = 0; // SNPs
		word *caseData = dataset->GetVarCase(OIDX, idx);
		word *ctrlData = dataset->GetVarCtrl(OIDX, idx);

		WordByte wb;
		for (uint32 i = 0; i < dataset->numWordCase; i++)
		{
			wb.w = caseData[i];
			contingencyCase[wb.b[0]]++;
			contingencyCase[wb.b[1]]++;
			contingencyCase[wb.b[2]]++;
			contingencyCase[wb.b[3]]++;
			contingencyCase[wb.b[4]]++;
			contingencyCase[wb.b[5]]++;
			contingencyCase[wb.b[6]]++;
			contingencyCase[wb.b[7]]++;
		}

		for (uint32 i = 0; i < dataset->numWordCtrl; i++)
		{
			wb.w = ctrlData[i];
			contingencyCtrl[wb.b[0]]++;
			contingencyCtrl[wb.b[1]]++;
			contingencyCtrl[wb.b[2]]++;
			contingencyCtrl[wb.b[3]]++;
			contingencyCtrl[wb.b[4]]++;
			contingencyCtrl[wb.b[5]]++;
			contingencyCtrl[wb.b[6]]++;
			contingencyCtrl[wb.b[7]]++;
		}
	}

	void OR_2x(varIdx idx)
	{
		const uint32 OIDX = 1; // Pair
		word *caseData = dataset->GetVarCase(OIDX, idx);
		word *ctrlData = dataset->GetVarCtrl(OIDX, idx);

		WordByte wb;
		for (uint32 i = 0; i < dataset->numWordCase; i++)
		{
			wb.w = epiCaseWord[OIDX - 1][i] | caseData[i];
			contingencyCase[wb.b[0]]++;
			contingencyCase[wb.b[1]]++;
			contingencyCase[wb.b[2]]++;
			contingencyCase[wb.b[3]]++;
			contingencyCase[wb.b[4]]++;
			contingencyCase[wb.b[5]]++;
			contingencyCase[wb.b[6]]++;
			contingencyCase[wb.b[7]]++;
		}

		for (uint32 i = 0; i < dataset->numWordCtrl; i++)
		{
			wb.w = epiCtrlWord[OIDX - 1][i] | ctrlData[i];
			contingencyCtrl[wb.b[0]]++;
			contingencyCtrl[wb.b[1]]++;
			contingencyCtrl[wb.b[2]]++;
			contingencyCtrl[wb.b[3]]++;
			contingencyCtrl[wb.b[4]]++;
			contingencyCtrl[wb.b[5]]++;
			contingencyCtrl[wb.b[6]]++;
			contingencyCtrl[wb.b[7]]++;
		}
	}

	void OR_3x(varIdx idx)
	{
		const uint32 OIDX = 2; // Triplet
		word *caseData = dataset->GetVarCase(OIDX, idx);
		word *ctrlData = dataset->GetVarCtrl(OIDX, idx);

		WordByte wb;
		for (uint32 i = 0; i < dataset->numWordCase; i++)
		{
			wb.w = epiCaseWord[OIDX - 1][i] | caseData[i];
			contingencyCase[wb.b[0]]++;
			contingencyCase[wb.b[1]]++;
			contingencyCase[wb.b[2]]++;
			contingencyCase[wb.b[3]]++;
			contingencyCase[wb.b[4]]++;
			contingencyCase[wb.b[5]]++;
			contingencyCase[wb.b[6]]++;
			contingencyCase[wb.b[7]]++;
		}

		for (uint32 i = 0; i < dataset->numWordCtrl; i++)
		{
			wb.w = epiCtrlWord[OIDX - 1][i] | ctrlData[i];
			contingencyCtrl[wb.b[0]]++;
			contingencyCtrl[wb.b[1]]++;
			contingencyCtrl[wb.b[2]]++;
			contingencyCtrl[wb.b[3]]++;
			contingencyCtrl[wb.b[4]]++;
			contingencyCtrl[wb.b[5]]++;
			contingencyCtrl[wb.b[6]]++;
			contingencyCtrl[wb.b[7]]++;
		}
	}

	void OR_4x(varIdx idx)
	{
		const uint32 OIDX = 3; // Quadlet
		word *caseData = dataset->GetVarCase(OIDX, idx);
		word *ctrlData = dataset->GetVarCtrl(OIDX, idx);

		WordByte wb;
		for (uint32 i = 0; i < dataset->numWordCase; i++)
		{
			wb.w = epiCaseWord[OIDX - 1][i] | caseData[i];
			contingencyCase[wb.b[0]]++;
			contingencyCase[wb.b[1]]++;
			contingencyCase[wb.b[2]]++;
			contingencyCase[wb.b[3]]++;
			contingencyCase[wb.b[4]]++;
			contingencyCase[wb.b[5]]++;
			contingencyCase[wb.b[6]]++;
			contingencyCase[wb.b[7]]++;
		}

		for (uint32 i = 0; i < dataset->numWordCtrl; i++)
		{
			wb.w = epiCtrlWord[OIDX - 1][i] | ctrlData[i];
			contingencyCtrl[wb.b[0]]++;
			contingencyCtrl[wb.b[1]]++;
			contingencyCtrl[wb.b[2]]++;
			contingencyCtrl[wb.b[3]]++;
			contingencyCtrl[wb.b[4]]++;
			contingencyCtrl[wb.b[5]]++;
			contingencyCtrl[wb.b[6]]++;
			contingencyCtrl[wb.b[7]]++;
		}
	}
	
	void ResetContigencyTable_1()
	{
		memset(contingencyCtrl, 0, 4 * sizeof(sampleIdx));
		memset(contingencyCase, 0, 4 * sizeof(sampleIdx));
	}

	void ResetContigencyTable_2()
	{
		memset(contingencyCtrl, 0, 16 * sizeof(sampleIdx));
		memset(contingencyCase, 0, 16 * sizeof(sampleIdx));
	}

	void ResetContigencyTable_3()
	{
		memset(contingencyCtrl, 0, 64 * sizeof(sampleIdx));
		memset(contingencyCase, 0, 64 * sizeof(sampleIdx));
	}

	void ResetContigencyTable_4()
	{
		memset(contingencyCtrl, 0, 256 * sizeof(sampleIdx));
		memset(contingencyCase, 0, 256 * sizeof(sampleIdx));
	}

	double Gini_1()
	{
		const uint32 entry = 3;
		double purity = 0;

		for (uint32 i = 0; i < entry; i++)
		{
			uint32 index = cti[i];
			double nCase = (double)contingencyCase[index];
			double nCtrl = (double)contingencyCtrl[index];
			double sum = nCase + nCtrl;
			if (sum)
				purity += (P2(nCase) + P2(nCtrl)) / (sum * dataset->numSample);
		}
		return purity;
	}

	double Gini_2()
	{
		const uint32 entry = 9;
		double purity = 0;

		for (uint32 i = 0; i < entry; i++)
		{
			uint32 index = cti[i];
			double nCase = (double)contingencyCase[index];
			double nCtrl = (double)contingencyCtrl[index];
			double sum = nCase + nCtrl;
			if (sum)
				purity += (P2(nCase) + P2(nCtrl)) / (sum * dataset->numSample);
		}
		return purity;
	}
	
	double Gini_3()
	{
		const uint32 entry = 27;
		double purity = 0;

		for (uint32 i = 0; i < entry; i++)
		{
			uint32 index = cti[i];
			double nCase = (double)contingencyCase[index];
			double nCtrl = (double)contingencyCtrl[index];
			double sum = nCase + nCtrl;
			if (sum)
				purity += (P2(nCase) + P2(nCtrl)) / (sum * dataset->numSample);
		}
		return purity;
	}

	double Gini_4()
	{
		const uint32 entry = 81;
		double purity = 0;

		for (uint32 i = 0; i < entry; i++)
		{
			uint32 index = cti[i];
			double nCase = (double)contingencyCase[index];
			double nCtrl = (double)contingencyCtrl[index];
			double sum = nCase + nCtrl;
			if (sum)
				purity += (P2(nCase) + P2(nCtrl)) / (sum * dataset->numSample);
		}
		return purity;
	}

	void Epi_1(uint32 id)
	{
		const uint32 OIDX = 0; // SNP
		threadIdx = id;

		AllocateThreadMemory();

		printf("Thread %4u starting ...\n", threadIdx);

		varIdx idx[1];

		for (idx[0] = 0; idx[0] < dataset->numVariable; idx[0]++)
		{
			uint32 pt = (idx[0] % args.numThreads);
			if (pt == threadIdx)
			{
#ifdef PTEST
				clock_t xc1 = clock();
#endif
				ResetContigencyTable_1();
#ifdef PTEST
				clock_t xc2 = clock();
#endif
				OR_1x(idx[0]);
#ifdef PTEST
				clock_t xc3 = clock();
#endif
				// compute purity
				double p = Gini_1();
#ifdef PTEST
				clock_t xc4 = clock();
#endif
				// report SNP combination if purity meet threshold
				if (args.printP[OIDX])
					if (p >= args.p[OIDX])
						fprintf(topPfile[threadIdx], "%f,%s\n", p, dataset->nameVariable[idx[0]]);

				// Save Purity to compute IG of next order
				if (args.saveP[OIDX])
					SnpPurity[idx[0]] = p;

				// compute Information Gained
				if (args.computeIG[OIDX])
				{
					double max_p = dataset->setPurity;

					double ig = p - max_p;

					// report SNP combination if IG meet threshold
					if (args.printIG[OIDX])
						if (ig >= args.ig[OIDX])
							fprintf(topIGfile[threadIdx], "%f,%s\n", ig, dataset->nameVariable[idx[0]]);

					// compute the best IG
					if (args.bestIG)
						dataset->results[threadIdx].Max_1(ig, p, idx);
				}
#ifdef PTEST
				clock_t xc5 = clock();
				elapse[1] += xc2 - xc1;
				elapse[2] += xc3 - xc2;
				elapse[3] += xc4 - xc3;
				elapse[4] += xc5 - xc4;
#endif
			}
		}

		printf("Thread %4u Finish\n", threadIdx);
		FreeThreadMemory();
	}
	
	void Epi_2(uint32 id)
	{
		const uint32 OIDX = 1; // Pair
		threadIdx = id;

		AllocateThreadMemory();

		printf("Thread %4u starting ...\n", threadIdx);

		varIdx idx[2];

		for (idx[0] = 0; idx[0] < (dataset->numVariable - OIDX); idx[0]++)
		{
			uint32 pt = (idx[0] % args.numThreads);
			if (pt == threadIdx)
			{
				OR_1(idx[0]);
				for (idx[1] = idx[0] + 1; idx[1] < dataset->numVariable; idx[1]++)
				{
#ifdef PTEST
					clock_t xc1 = clock();
#endif
					ResetContigencyTable_2();
#ifdef PTEST
					clock_t xc2 = clock();
#endif
					OR_2x(idx[1]);
#ifdef PTEST
					clock_t xc3 = clock();
#endif
					// compute purity
					double p = Gini_2();
#ifdef PTEST
					clock_t xc4 = clock();
#endif
					// report SNP combination if purity meet threshold
					if (args.printP[OIDX])
						if (p >= args.p[OIDX])
							fprintf(topPfile[threadIdx], "%f,%s,%s\n", p, dataset->nameVariable[idx[0]], dataset->nameVariable[idx[1]]);

					// Save Purity to compute IG of next order
					if (args.saveP[OIDX])
						PairPurity[idx[0]][idx[1]] = p;

					// compute Information Gained
					if (args.computeIG[OIDX])
					{
						double max_p = (SnpPurity[idx[1]] > SnpPurity[idx[0]]) ? SnpPurity[idx[1]] : SnpPurity[idx[0]];
		
						double ig = p - max_p;

						// report SNP combination if IG meet threshold
						if (args.printIG[OIDX])
							if (ig >= args.ig[OIDX])
								fprintf(topIGfile[threadIdx], "%f,%s,%s\n", ig, dataset->nameVariable[idx[0]], dataset->nameVariable[idx[1]]);

						// compute the best IG
						if (args.bestIG)
							dataset->results[threadIdx].Max_2(ig, p, idx);
					}
#ifdef PTEST
					clock_t xc5 = clock();
					elapse[1] += xc2 - xc1;
					elapse[2] += xc3 - xc2;
					elapse[3] += xc4 - xc3;
					elapse[4] += xc5 - xc4;
#endif
				}
			}
		}
		printf("Thread %4u Finish\n", threadIdx);
		FreeThreadMemory();
	}

	void Epi_3(uint32 id)
	{
		const uint32 OIDX = 2; // Triplet
		threadIdx = id;

		AllocateThreadMemory();

		printf("Thread %4u starting ...\n", threadIdx);

		varIdx idx[3];

		for (idx[0] = 0; idx[0] < (dataset->numVariable - OIDX); idx[0]++)
		{
			uint32 pt = (idx[0] % args.numThreads);
			if (pt == threadIdx)
			{
				OR_1(idx[0]);
				for (idx[1] = idx[0] + 1; idx[1] < (dataset->numVariable - (OIDX - 1)); idx[1]++)
				{
					OR_2(idx[1]);
					for (idx[2] = idx[1] + 1; idx[2] < dataset->numVariable; idx[2]++)
					{
#ifdef PTEST
						clock_t xc1 = clock();
#endif
						ResetContigencyTable_3();
#ifdef PTEST
						clock_t xc2 = clock();
#endif
						OR_3x(idx[2]);
#ifdef PTEST
						clock_t xc3 = clock();
#endif
						// compute purity
						double p = Gini_3();
#ifdef PTEST
						clock_t xc4 = clock();
#endif
						// report SNP combination if purity meet threshold
						if (args.printP[OIDX])
							if (p >= args.p[OIDX])
								fprintf(topPfile[threadIdx], "%f,%s,%s,%s\n", p, dataset->nameVariable[idx[0]], dataset->nameVariable[idx[1]], dataset->nameVariable[idx[2]]);

						// Save Purity to compute IG of next order
						if (args.saveP[OIDX])
							tripletPurity[idx[0]][idx[1]][idx[2]] = p;

						// compute Information Gained
						if (args.computeIG[OIDX])
						{
							double max_p = (PairPurity[idx[0]][idx[1]] > PairPurity[idx[0]][idx[2]]) ? PairPurity[idx[0]][idx[1]] : PairPurity[idx[0]][idx[2]];
							max_p = (PairPurity[idx[1]][idx[2]] > max_p) ? PairPurity[idx[1]][idx[2]] : max_p;

							double ig = p - max_p;

							// report SNP combination if IG meet threshold
							if (args.printIG[OIDX])
								if (ig >= args.ig[OIDX])
									fprintf(topIGfile[threadIdx], "%f,%s,%s,%s\n", ig, dataset->nameVariable[idx[0]], dataset->nameVariable[idx[1]], dataset->nameVariable[idx[2]]);

							// compute the best IG
							if (args.bestIG)
								dataset->results[threadIdx].Max_3(ig, p, idx);
						}
#ifdef PTEST
						clock_t xc5 = clock();
						elapse[1] += xc2 - xc1;
						elapse[2] += xc3 - xc2;
						elapse[3] += xc4 - xc3;
						elapse[4] += xc5 - xc4;
#endif
					}
				}
			}
		}
		printf("Thread %4u Finish\n", threadIdx);
		FreeThreadMemory();
	}

	void Epi_4(uint32 id)
	{
		const uint32 OIDX = 3; // Quadlet
		threadIdx = id;

		AllocateThreadMemory();

		printf("Thread %4u starting ...\n", threadIdx);

		varIdx idx[4];

		for (idx[0] = 0; idx[0] < (dataset->numVariable - OIDX); idx[0]++)
		{
			uint32 pt = (idx[0] % args.numThreads);
			if (pt == threadIdx)
			{
				OR_1(idx[0]);
				for (idx[1] = idx[0] + 1; idx[1] < (dataset->numVariable - (OIDX - 1)); idx[1]++)
				{
					OR_2(idx[1]);
					for (idx[2] = idx[1] + 1; idx[2] < (dataset->numVariable - (OIDX - 2)); idx[2]++)
					{
						OR_3(idx[2]);
						for (idx[3] = idx[2] + 1; idx[3] < dataset->numVariable; idx[3]++)
						{
#ifdef PTEST
							clock_t xc1 = clock();
#endif
							ResetContigencyTable_4();
#ifdef PTEST
							clock_t xc2 = clock();
#endif
							OR_4x(idx[3]);
#ifdef PTEST
							clock_t xc3 = clock();
#endif
							// compute purity
							double p = Gini_4();
#ifdef PTEST
							clock_t xc4 = clock();
#endif

							// report SNP combination if purity meet threshold
							if (args.printP[OIDX])
								if (p >= args.p[OIDX])
									fprintf(topPfile[threadIdx], "%f,%s,%s,%s,%s\n", p, dataset->nameVariable[idx[0]], dataset->nameVariable[idx[1]], dataset->nameVariable[idx[2]], dataset->nameVariable[idx[3]]);

							// compute Information Gained
							if (args.computeIG[OIDX])
							{
								double max_p = (tripletPurity[idx[0]][idx[1]][idx[2]] > tripletPurity[idx[0]][idx[1]][idx[3]]) ? tripletPurity[idx[0]][idx[1]][idx[2]] : tripletPurity[idx[0]][idx[1]][idx[3]];
								max_p = (tripletPurity[idx[0]][idx[2]][idx[3]] > max_p) ? tripletPurity[idx[0]][idx[2]][idx[3]] : max_p;
								max_p = (tripletPurity[idx[1]][idx[2]][idx[3]] > max_p) ? tripletPurity[idx[1]][idx[2]][idx[3]] : max_p;

								double ig = p - max_p;

								// report SNP combination if IG meet threshold
								if (args.printIG[OIDX])
									if (ig >= args.ig[OIDX])
										fprintf(topIGfile[threadIdx], "%f,%s,%s,%s,%s\n", ig, dataset->nameVariable[idx[0]], dataset->nameVariable[idx[1]], dataset->nameVariable[idx[2]], dataset->nameVariable[idx[3]]);

								// compute the best IG
								if (args.bestIG)
									dataset->results[threadIdx].Max_4(ig, p, idx);
							}
#ifdef PTEST
							clock_t xc5 = clock();
							elapse[1] += xc2 - xc1;
							elapse[2] += xc3 - xc2;
							elapse[3] += xc4 - xc3;
							elapse[4] += xc5 - xc4;
#endif
						}
					}
				}
			}
		}

		printf("Thread %4u Finish\n", threadIdx);
		FreeThreadMemory();
	}

	void MultiThread(void *(*threadFunction) (void *))
	{
		pthread_t *threads = new pthread_t[args.numThreads];
		ThreadData *td = new ThreadData[args.numThreads];
		for (uint32 i = 0; i < args.numThreads; i++)
		{
			td[i].epiStat = (void *) new EpiStat(this);
			td[i].id = i;
		}
		for (uint32 i = 0; i < args.numThreads; i++)
		{
			pthread_create(&threads[i], NULL, threadFunction, &td[i]);
		}
		for (uint32 i = 0; i < args.numThreads; i++)
		{
			pthread_join(threads[i], NULL);
		}

		
	}

	void Run()
	{
		for (int i = 0; i < MAX_ORDER; i++)
		{
			if (args.computeP[i])
			{
				time_t begin = time(NULL);
				printf("\n\n>>>>>>>>>> Process %u-SNP combinations\n", i + 1);

				OpenFiles(i);
				MultiThread(threadFunction[i]);
				CloseFiles(i);

				time_t end = time(NULL);
				double time_spent = difftime(end, begin);
				printf("\n\n<<<<<<<<< Prosess %u-SNP combinations takes %10.0f seconds\n\n", i+1, time_spent);
			}
		}
		// merge thread files
		#ifndef _MSC_VER
		{
			char* cmd = new char[1024];
			NULL_CHECK(cmd);
			for (uint32 order = 0; order < MAX_ORDER; order++)
			{
				char header[1024];
				switch (order)
				{
				case 0:	sprintf(header, "SNP_A"); break;
				case 1:	sprintf(header, "SNP_A,SNP_B"); break;
				case 2:	sprintf(header, "SNP_A,SNP_B,SNP_C"); break;
				case 3:	sprintf(header, "SNP_A,SNP_B,SNP_C,SNP_D"); break;
				}

				char sortCmd[100];
				if (args.sort)
					sprintf(sortCmd, "%s", "| sort -g -r -k1,1 -t ','");
				else
					sprintf(sortCmd, " ");

				if (args.printP[order])
				{
					// create a merged output file
					sprintf(cmd, "cat %s.Purity.%u.*.csv %s | awk 'BEGIN{print(\"Purity,%s\")}{print}' > %s.Purity.%u.csv", args.output, order, sortCmd, header, args.output, order);
					printf("\n>>> %s\n", cmd);
					if (system(cmd) == -1)
						ERROR("Cannot merge output files");
					sprintf(cmd, "rm %s.Purity.%u.*.csv", args.output, order);
					printf("\n>>> %s\n", cmd);
					if (system(cmd) == -1)
						ERROR("Cannot delete temp files");
				}
				if (args.printIG[order])
				{
					// create a merged output file
					sprintf(cmd, "cat %s.IG.%u.*.csv %s | awk 'BEGIN{print(\"IG,%s\")}{print}' > %s.IG.%u.csv", args.output, order, sortCmd, header, args.output, order);
					printf("\n>>> %s\n", cmd);
					if (system(cmd) == -1)
						ERROR("Cannot merge output files");
					sprintf(cmd, "rm %s.IG.%u.*.csv", args.output, order);
					printf("\n>>> %s\n", cmd);
					if (system(cmd) == -1)
						ERROR("Cannot delete temp files");
				}
			}
			delete[]cmd;
		}
		#endif
		if (args.bestIG)
		{
			printf("\n\nAggregate result from threads\n");
			for (uint32 i = 1; i < args.numThreads; i++)
				dataset->results[0].Max(dataset->results[i]);

			char* fn = new char[strlen(args.output) + 20];
			NULL_CHECK(fn);
			sprintf(fn, "%s.bestIG.csv", args.output);
			dataset->results->toCSV(fn, dataset->nameVariable);
			delete[]fn;
		}

		FreePurity(dataset->numVariable, args);
		dataset->FreeMemory(args);
	}
};

void *EpiThread_1(void *t)
{
	ThreadData *td = (ThreadData *)t;
	EpiStat *epiStat = (EpiStat *)td->epiStat;
	epiStat->Epi_1(td->id);
	return NULL;
}

void *EpiThread_2(void *t)
{
	ThreadData *td = (ThreadData *)t;
	EpiStat *epiStat = (EpiStat *)td->epiStat;
	epiStat->Epi_2(td->id);
	return NULL;
}

void *EpiThread_3(void *t)
{
	ThreadData *td = (ThreadData *)t;
	EpiStat *epiStat = (EpiStat *)td->epiStat;
	epiStat->Epi_3(td->id);
	return NULL;
}

void *EpiThread_4(void *t)
{
	ThreadData *td = (ThreadData *)t;
	EpiStat *epiStat = (EpiStat *)td->epiStat;
	epiStat->Epi_4(td->id);
	return NULL;
}

int main(int argc, char *argv[])
{
	#ifdef PTEST
		for (uint32 i = 0; i < 100; i++)
			elapse[i] = 0;
		clock_t xc1 = clock();
	#endif

	ARGS args;
	args.Parse(argc, argv);
	//args.Print();
	
	Dataset dataset;
	dataset.ReadDataset(args.input);
	dataset.Init(args);

	AllocatePurity(dataset.numVariable, args);

	EpiStat epiStat;
	epiStat.Init(&dataset, args, EpiThread_1, EpiThread_2, EpiThread_3, EpiThread_4);

	epiStat.Run();

	printf("\n=============Finish=============\n\n\n");

	#ifdef PTEST
		clock_t xc2 = clock();
		elapse[0] = xc2 - xc1;
		clock_t other = elapse[0];
		for (uint32 i = 1; i < 100; i++)
		{
			elapse[i] /= args.numThreads;
			other -= elapse[i];
		}
		
		printf("\n Total: %15u", elapse[0]);
		printf("\n Clear: %15u \t %5.2f %%", elapse[1], ((double)elapse[1] *100) / elapse[0]);
		printf("\n Count: %15u \t %5.2f %%", elapse[2], ((double)elapse[2] *100) / elapse[0]);
		printf("\n Gini : %15u \t %5.2f %%", elapse[3], ((double)elapse[3] *100) / elapse[0]);
		printf("\n Rest : %15u \t %5.2f %%", elapse[4], ((double)elapse[4] *100) / elapse[0]);
		printf("\n Other: %15u \t %5.2f %%", other,     ((double)other     *100) / elapse[0]);
		printf("\n\n");
	#endif
	return 0;
}

