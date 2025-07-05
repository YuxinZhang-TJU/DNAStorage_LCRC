/************************************
 * Copyright (C) 2025 Tianjin University (TJU)
 *
 * @brief: DNA encoder for large DNA data storage
 *
 *************************************/
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <cmath>
#include <cassert>
#include <cstdint>

#include "compositeCode.h"
#include "utils.h"

#define COMPONENT_N 5 // Number of component codes

#define LDPC_N 22680
#define LDPC_K 7560

const char *src_encoded_data_file = "srcdat/Tang_Poem_codeword.txt";
const unsigned ldpc_num = 1;

int main(int argc, char *argv[])
{
    setvbuf(stdout, NULL, _IONBF, 0); // none buffer stdout

    char filename_largeDNA[500] = {0};

    unsigned cLen[COMPONENT_N] = {31, 35, 43, 47, 59};
    int cPhase_init[COMPONENT_N] = {0};

    // Lengths of component codes
    unsigned start_phase = 156802;
    for (int i = 0; i < COMPONENT_N; i++)
        *(cPhase_init + i) = start_phase % *(cLen + i);

    if (argc != 2)
    {
        fprintf(stderr, "Usage: "
                        "%s [output DNA fasta]\n\n",
                argv[0]);
        return 1;
    }
    sscanf(argv[1], "%s", filename_largeDNA);

    time_t start_time, end_time;
    char timebuff[100];
    time(&start_time);
    strftime(timebuff, 100, "%Y-%m-%d %H:%M:%S ", localtime(&start_time));
    printf("\nStart running: %s\n\n", timebuff);

    if (!isCopirme(cLen, COMPONENT_N))
    {
        fprintf(stderr, "\nPeriods of components must be coprime!\n\n");
        return 1;
    }

    unsigned period = 1;
    for (int i = 0; i < COMPONENT_N; i++)
        period *= *(cLen + i);

    unsigned largeDNA_len = ldpc_num * LDPC_N;

    printf("Component code num = %d, periods: ", COMPONENT_N);
    for (int i = 0; i < COMPONENT_N; i++)
        printf("%d ", *(cLen + i));
    printf("\nFull length composite code period: %u\n", period);
    printf("LDPC(%d, %d) num = %u, large DNA len = %u\n\n", LDPC_N, LDPC_K, ldpc_num, largeDNA_len);

    printf("----Step 1: Generation of composite code----\n\n");
    //========================== Generation of composite code =======================//
    int **c = (int **)malloc(COMPONENT_N * sizeof(int *));
    for (int i = 0; i < COMPONENT_N; i++)
    {
        *(c + i) = (int *)alloc_init(*(cLen + i), sizeof(int));

        if (PNSequenceGenerator(*(c + i), *(cLen + i), 0))
            return 1;
    }

    char *compositeCode = (char *)alloc_init(largeDNA_len, sizeof(char));
    if (createCompositeCode_majority(compositeCode, largeDNA_len, COMPONENT_N, cPhase_init, cLen, c))
        return 1;
    printf("Length of composite code used: %d\n\n", largeDNA_len);

    printf("----Step 2: LDPC encoding and DNA transcoding----\n\n");

    //========================== LDPC(22680, 7560) encoding =======================//
    char src_codeword[LDPC_N] = {0};
    char codeword_permutated[LDPC_N] = {0};
    FILE *fpr_src = NULL;
    if ((fpr_src = fopen(src_encoded_data_file, "r")) == NULL)
    {
        fprintf(stderr, "Failed to open %s\n", src_encoded_data_file);
        exit(1);
    }
    fscanf(fpr_src, "%s\n", src_codeword);
    fclose(fpr_src);

    unsigned *randVal = (unsigned *)alloc_init(32767, sizeof(unsigned)); // 2^15
    unsigned *permutation = (unsigned *)alloc_init(LDPC_N, sizeof(unsigned));
    PRNG(randVal, 15, 32767, 64);
    for (unsigned i = 0, nn = 0; i < 32767; ++i)
    {
        if (*(randVal + i) <= LDPC_N)
        {
            *(permutation + nn) = *(randVal + i) - 1;
            ++nn;
        }
    }

    char *encoded_largeDNA = (char *)alloc_init(largeDNA_len, sizeof(char));

    //========================== Transcoding into DNA =======================//
    for (unsigned tb = 0; tb < LDPC_N; ++tb)
    {
        *(codeword_permutated + tb) = *(src_codeword + *(permutation + tb));

        if (*(codeword_permutated + tb) == '0' &&
            *(compositeCode + tb) == 1)
        {
            *(encoded_largeDNA + tb) = 'A';
        }
        else if (*(codeword_permutated + tb) == '0' &&
                 *(compositeCode + tb) == -1)
        {
            *(encoded_largeDNA + tb) = 'T';
        }
        else if (*(codeword_permutated + tb) == '1' &&
                 *(compositeCode + tb) == 1)
        {
            *(encoded_largeDNA + tb) = 'G';
        }
        else if (*(codeword_permutated + tb) == '1' &&
                 *(compositeCode + tb) == -1)
        {
            *(encoded_largeDNA + tb) = 'C';
        }
    }

    FILE *fpw_fa = NULL;
    if ((fpw_fa = fopen(filename_largeDNA, "w")) == NULL)
    {
        fprintf(stderr, "Failed to open %s!\n", filename_largeDNA);
        exit(1);
    }
    fprintf(fpw_fa, ">1\n");
    fwrite(encoded_largeDNA, 1, largeDNA_len, fpw_fa);
    fprintf(fpw_fa, "\n");
    fclose(fpw_fa);

    printf("Encoding finished!\n"
           "The resulting large DNA (%u bp) is written to %s\n",
           largeDNA_len, filename_largeDNA);

    time(&end_time);
    strftime(timebuff, 100, "%Y-%m-%d %H:%M:%S ", localtime(&end_time));
    printf("\nEnd running: %s\n", timebuff);
    printf("\nElapsed time: %.2lf sec\n\n", difftime(end_time, start_time));
    return 0;
}
