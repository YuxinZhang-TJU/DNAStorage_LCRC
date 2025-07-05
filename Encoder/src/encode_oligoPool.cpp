/************************************
 * Copyright (C) 2025 Tianjin University (TJU)
 *
 * @brief: DNA encoder for oligo pool data storage
 *
 *************************************/
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <cmath>
#include <cassert>
#include <cstdint>

#include "ldpc_encode.h"
#include "compositeCode.h"
#include "utils.h"

#define COMPONENT_N 5 // Number of component codes

#define RAW_DATA_BYTE 195750

#define OLIGO_NUM 11745
#define OLIGO_LEN 200
#define PAYLOAD_LEN 160
#define PRIMER_LEN 20

#define LDPC_K 54000
#define LDPC_N 64800
#define LDPC_M 10800
#define CHECK_DEGREE 20

#define MASK_LEN 10800 * 20
#define BLOCKCODE_NUM 29

#define MSEQ1 65535
#define MSEQ2 31

const char *primerL = "ATAATTGGCTCCTGCTTGCA";
const char *primerR = "AATGTAGGCGGAAAGTGCAA";

const char *img1 = "srcdat/Grand_Canyon.jpg";
const char *img2 = "srcdat/Mount_Everest.jpg";

const char *enc_config_file = "config/ldpc_enc";

void ascii2bit(char *outBit, char symbol);
int read_file_and_write_bit(char *bitstream);

int main(int argc, char *argv[])
{
    setvbuf(stdout, NULL, _IONBF, 0);

    char filename_oligos[500] = {0};

    int primer_len = PRIMER_LEN;

    // Lengths of component codes
    unsigned cLen[COMPONENT_N] = {31, 35, 43, 47, 59};
    int cPhase_init[COMPONENT_N] = {0};

    if (argc != 2)
    {
        fprintf(stderr, "Usage: "
                        "%s [output DNA fasta]\n\n",
                argv[0]);
        return 1;
    }
    sscanf(argv[1], "%s", filename_oligos);

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

    unsigned mSeq1Len = MSEQ1;
    unsigned mSeq2Len = MSEQ2;
    char mSeq1[MSEQ1];
    char mSeq2[MSEQ2];
    FILE *fpr = NULL;
    if ((fpr = fopen("config/m65535.txt", "rb")) == NULL)
    {
        fprintf(stderr, "Could not open file config/m65535.txt!\n");
        return 1;
    }
    fread(mSeq1, sizeof(char), mSeq1Len, fpr);
    fclose(fpr);
    for (int i = 0; i < mSeq1Len; i++)
    {
        *(mSeq1 + i) -= '0';
    }
    if ((fpr = fopen("config/m31.txt", "rb")) == NULL)
    {
        fprintf(stderr, "Could not open file config/m31.txt!\n");
        return 1;
    }
    fread(mSeq2, sizeof(char), mSeq2Len, fpr);
    fclose(fpr);
    for (int i = 0; i < mSeq2Len; i++)
    {
        *(mSeq2 + i) -= '0';
    }

    int *encoder_mask = (int *)calloc(LDPC_M * CHECK_DEGREE, sizeof(int));
    FILE *fpr_config = NULL;
    if ((fpr_config = fopen(enc_config_file, "rb")) == NULL)
    {
        fprintf(stderr, "Could not open config file '%s'!\n", enc_config_file);
        exit(1);
    }
    fread(encoder_mask, 4, LDPC_M * CHECK_DEGREE, fpr_config);
    fclose(fpr_config);

    unsigned period = 1;
    for (int i = 0; i < COMPONENT_N; i++)
        period *= *(cLen + i);
    assert(LDPC_N * BLOCKCODE_NUM <= period);

    unsigned compositeCode_used_len = BLOCKCODE_NUM * LDPC_N;

    printf("Component code num = %d, periods: ", COMPONENT_N);
    for (int i = 0; i < COMPONENT_N; i++)
        printf("%d ", *(cLen + i));
    printf("\nFull length composite code period: %u\n", period);
    printf("LDPC(%d, %d) num = %u, encoded data len = %u\n\n", LDPC_N, LDPC_K, BLOCKCODE_NUM, compositeCode_used_len);

    char cblk[LDPC_M];

    char *rawData = (char *)alloc_init(BLOCKCODE_NUM * LDPC_K, sizeof(char));
    char *rawData_ptr = NULL;

    char **codeword = (char **)malloc(BLOCKCODE_NUM * sizeof(char *));
    for (int i = 0; i < BLOCKCODE_NUM; ++i)
    {
        codeword[i] = (char *)alloc_init(LDPC_N, sizeof(char));
    }
    char *scrambleData = (char *)alloc_init(LDPC_K, sizeof(char));
    char *permutedCode = (char *)alloc_init(LDPC_N, sizeof(char));

    char *encData = (char *)alloc_init(LDPC_N * BLOCKCODE_NUM, sizeof(char));

    char **oligos = (char **)malloc(OLIGO_NUM * sizeof(char *));
    for (int i = 0; i < OLIGO_NUM; ++i)
        oligos[i] = (char *)alloc_init(OLIGO_LEN + 1, sizeof(char));

    read_file_and_write_bit(rawData);

    printf("----Step 1: Generation of composite code----\n\n");

    //========================== Generation of composite code =======================//
    int **c = (int **)malloc(COMPONENT_N * sizeof(int *)); // Component sequences
    for (int i = 0; i < COMPONENT_N; i++)
    {
        *(c + i) = (int *)alloc_init(*(cLen + i), sizeof(int));

        if (PNSequenceGenerator(*(c + i), *(cLen + i), 0))
            return 1;
    }

    char *compositeCode = (char *)alloc_init(compositeCode_used_len, sizeof(char));
    if (createCompositeCode_majority(compositeCode, compositeCode_used_len, COMPONENT_N, cPhase_init, cLen, c))
        return 1;
    char *compositeCode_bin = (char *)alloc_init(compositeCode_used_len, sizeof(char));
    for (unsigned ti = 0; ti < compositeCode_used_len; ++ti)
        *(compositeCode_bin + ti) = *(compositeCode + ti) == 1 ? 0 : 1;
    printf("Length of composite code used: %d\n\n", compositeCode_used_len);

    printf("----Step 2: LDPC encoding and DNA transcoding----\n\n");

    //========================== LDPC(64800, 54000) encoding =======================//
    unsigned *randVal = (unsigned *)alloc_init(65535, sizeof(unsigned));
    unsigned *permutation = (unsigned *)alloc_init(LDPC_N, sizeof(unsigned));
    PRNG(randVal, 16, 65535, 64);
    for (unsigned i = 0, nn = 0; i < 65535; ++i)
    {
        if (*(randVal + i) <= LDPC_N)
        {
            *(permutation + nn) = *(randVal + i) - 1;
            ++nn;
        }
    }

    for (int i = 0; i < BLOCKCODE_NUM; ++i)
    {
        rawData_ptr = rawData + LDPC_K * i;
        for (int j = 0; j < LDPC_K; ++j)
        {
            *(scrambleData + j) = (*(rawData_ptr + j) +
                                   *(mSeq1 + ((unsigned)LDPC_K * i + j) % mSeq1Len) +
                                   *(mSeq2 + ((unsigned)LDPC_K * i + j) % mSeq2Len)) %
                                  2;
        }

        ldpc_encode(scrambleData, cblk, encoder_mask);

        for (int j = 0; j < LDPC_M; ++j)
        {
            *(*(codeword + i) + j) = cblk[j];
        }
        for (int j = 0; j < LDPC_K; ++j)
        {
            *(*(codeword + i) + LDPC_M + j) = *(scrambleData + j);
        }

        for (int j = 0; j < LDPC_N; ++j)
            *(permutedCode + j) = *(*(codeword + i) + *(permutation + j));

        memcpy(encData + LDPC_N * i, permutedCode, sizeof(char) * LDPC_N);
    }

    //========================== Transcoding into DNA =======================//
    char *upper_payload = NULL;
    char *lower_composite = NULL;
    for (int i = 0; i < OLIGO_NUM; ++i)
    {
        upper_payload = encData + PAYLOAD_LEN * i;
        lower_composite = compositeCode_bin + PAYLOAD_LEN * i;

        for (int j = 0; j < PAYLOAD_LEN; ++j)
        {
            if (*(upper_payload + j) == 0 & *(lower_composite + j) == 0)
                *(*(oligos + i) + primer_len + j) = 'A';
            else if (*(upper_payload + j) == 0 & *(lower_composite + j) == 1)
                *(*(oligos + i) + primer_len + j) = 'T';
            else if (*(upper_payload + j) == 1 & *(lower_composite + j) == 0)
                *(*(oligos + i) + primer_len + j) = 'G';
            else if (*(upper_payload + j) == 1 & *(lower_composite + j) == 1)
                *(*(oligos + i) + primer_len + j) = 'C';
        }

        for (int j = 0; j < primer_len; ++j)
        {
            *(*(oligos + i) + j) = *(primerL + j);
            *(*(oligos + i) + primer_len + PAYLOAD_LEN + j) = *(primerR + j);
        }
    }

    FILE *fpw = NULL;
    if ((fpw = fopen(filename_oligos, "wb")) == NULL)
    {
        fprintf(stderr, "Failed to open %s\n", filename_oligos);
        exit(1);
    }
    for (int i = 0; i < OLIGO_NUM; ++i)
    {
        fprintf(fpw, ">%d\n", i + 1);
        fwrite(*(oligos + i), 1, OLIGO_LEN, fpw);
        fwrite("\n", 1, 1, fpw);
    }
    fclose(fpw);

    printf("Encoding finished!\n"
           "The resulting oligos (%d x %d nt) are written to %s\n",
           OLIGO_NUM, OLIGO_LEN, filename_oligos);

    time(&end_time);
    strftime(timebuff, 100, "%Y-%m-%d %H:%M:%S ", localtime(&end_time));
    printf("\nEnd running: %s\n", timebuff);
    printf("Elapsed time: %.2lf sec\n\n", difftime(end_time, start_time));

    return 0;
}

void ascii2bit(char *outBit, char symbol)
{
    for (int i = 0; i < 8; i++)
    {
        *(outBit + i) = ((symbol >> (7 - i)) & 1);
    }
}

int read_file_and_write_bit(char *bitstream)
{
    FILE *fpr = NULL;
    size_t totSize = 0;
    size_t totbits = 0;
    char *pbuff;

    if ((fpr = fopen(img1, "rb")) == NULL)
    {
        fprintf(stderr, "Faile to open %s\n", img1);
        exit(1);
    }
    fseek(fpr, 0, SEEK_END);
    long fileSize = ftell(fpr);
    totSize += fileSize;
    pbuff = (char *)calloc(fileSize, sizeof(char));
    rewind(fpr);
    fread(pbuff, sizeof(char), fileSize, fpr);
    fclose(fpr);
    fpr = NULL;
    printf("Data volume of %s: %ld bytes\n", img1, fileSize);

    char bitBuff[8] = {0};
    for (long i = 0; i < fileSize; ++i)
    {
        ascii2bit(bitBuff, *(pbuff + i));

        memcpy(bitstream + totbits, bitBuff, 8);
        totbits += 8;
    }
    free(pbuff);

    if ((fpr = fopen(img2, "rb")) == NULL)
    {
        fprintf(stderr, "Faile to open %s\n", img2);
        exit(1);
    }
    fseek(fpr, 0, SEEK_END);
    fileSize = ftell(fpr);
    totSize += fileSize;
    pbuff = (char *)calloc(fileSize, sizeof(char));
    rewind(fpr);
    fread(pbuff, sizeof(char), fileSize, fpr);
    fclose(fpr);
    fpr = NULL;

    for (long i = 0; i < fileSize; ++i)
    {
        ascii2bit(bitBuff, *(pbuff + i));

        memcpy(bitstream + totbits, bitBuff, 8);
        totbits += 8;
    }
    free(pbuff);

    printf("Data volume of %s: %ld bytes\n", img2, fileSize);
    printf("Total data volume: %ld bytes\n", totSize);
    printf("Padding with 0: %ld bytes\n\n", RAW_DATA_BYTE - totSize);

    assert(totSize <= RAW_DATA_BYTE);

    return 0;
}
