/* Date : 5th July, 2014 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Maximum length of FASTA sequence ID */
#define BUFFER 1024

/* FASTA type sequence data structure. */
typedef struct {
    char _id[BUFFER];
    /* Sequence length of size upto 10*1024 is supported. For more length
     * increase this size or use malloc/realloc
     */
    char _sequence[10*BUFFER];
} fasta_entry;

int calculate_count(char *val, FILE *out) {
    size_t len = strlen(val);
    size_t i;

    /* The four indexes store count of A,T,G and C respectively */
    int count[4] = {0,0,0,0};
    for(i = 0; i < len; ++i) {
        if(val[i] == 'A' || val[i] == 'a') count[0]++;
        if(val[i] == 'T' || val[i] == 't') count[1]++;
        if(val[i] == 'G' || val[i] == 'g') count[2]++;
        if(val[i] == 'C' || val[i] == 'c') count[3]++;
    }

    /* Write output to file. */
    fprintf(out, "Total number of residues :\n");
    fprintf(out, "A: %d,\nT: %d,\nG: %d,\nC: %d\n", count[0], count[1],\
            count[2], count[3]);

    return 0;
}
int main(int argc, char **argv) {
    /* Throw an error if filename is not supplied */
    /*
    if(argc!=2) {
        fprintf(stderr, "Incorrect Command. Please provide file.\n");
        exit(-1);
    }
    */

    /* I/O Variables and data strctures */
    FILE *in, *out;
    char line[BUFFER];
    int counter = 0;
    size_t i;
    fasta_entry *entries;
    
    /* Space for 1024 FASTA entries is allocated. Make necessary changes
     * before trying for higher number of sequences.*/
    entries = (fasta_entry *) malloc(sizeof(fasta_entry)*BUFFER);

    in = fopen("input.txt", "r");
    out = fopen("result.txt", "w");

    /* Throw error if cannot open input file. */
    if(in == NULL) {
        fprintf(stderr, "Error opening file.\n");
    }

    /* Reading input file line by line until NULL.
     * IMPORTANT: The code is written assuming file format to be UNIX 
     * compatible i.e, newlines are given by \n and not \r\n like 
     * windows. Make sure you run it on correct platform.
     */
    while(fgets(line, BUFFER, in) != NULL) {
        size_t ln = strlen(line);
        if(line[0] == '\n' && line[1] == '>') {
            strcpy(entries[counter]._id, line);
            ++counter;
        }
        else if(line[0] == '>') {
            strcpy(entries[counter]._id, line);
            ++counter;
        }
        else {
            if(line[ln-1] == '\n') line[ln-1] = '\0';
            strcat(entries[counter-1]._sequence, line);
        }
    }

    /* Give a successful reading hint at stdout. */
    fprintf(stdout, "Total %d entries recorded.\n", counter);

    /* Write to file. */
    for(i = 0; i < counter; ++i) {
        fprintf(out, "ID: %s", entries[i]._id);
        calculate_count(entries[i]._sequence, out);
        fprintf(out, "\n");
    }

    /* Close all open file handles, free memory. */
    fclose(in);
    fclose(out);
    free(entries);

    return 0;
}
