#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define MAX_LINE 256
#define ENZYME_COUNT 60

typedef struct {
    char name[16];
    char pattern[16];
} Enzyme;

Enzyme enzymes[ENZYME_COUNT] = {
    {"EcoRI", "GAATTC"}, {"BamHI", "GGATCC"}, {"HindIII","AAGCTT"},
    {"TaqI", "TCGA"}, {"HaeIII", "GGCC"}, {"AluI", "AGCT"},
    {"AccI", "GTMKAC"}, {"AflII", "CTTAAG"}, {"ApaI", "GGGCCC"},
    {"BglII", "AGATCT"}, {"ClaI", "ATCGAT"}, {"DraI", "TTTAAA"},
    {"EcoRV", "GATATC"}, {"FokI", "GGATG"}, {"HhaI", "GCGC"},
    {"HincII", "GTYRAC"}, {"KpnI", "GGTACC"}, {"MboI", "GATC"},
    {"MspI", "CCGG"}, {"NcoI", "CCATGG"}, {"NdeI", "CATATG"},
    {"NotI", "GCGGCCGC"}, {"NsiI", "ATGCAT"}, {"PstI", "CTGCAG"},
    {"PvuII", "CAGCTG"}, {"SacI", "GAGCTC"}, {"SalI", "GTCGAC"},
    {"ScaI", "AGTACT"}, {"SmaI", "CCCGGG"}, {"SpeI", "ACTAGT"},
    {"SphI", "GCATGC"}, {"SspI", "AATATT"}, {"StuI", "AGGCCT"},
    {"Tth111I", "GACNNNGTC"}, {"XbaI", "TCTAGA"}, {"XhoI", "CTCGAG"},
    {"XmaI", "CCCGGG"}, {"BsaI", "GGTCTC"}, {"BsmBI", "CGTCTC"},
    {"BsrGI", "TGTACA"}, {"BstEII", "GGTNACC"}, {"EagI", "CGGCCG"},
    {"AatII", "GACGTC"}, {"BbvCI", "CCTCAGC"}, {"BsiWI", "CGTACG"},
    {"BspEI", "TCCGGA"}, {"BsrI", "ACTGG"}, {"BstAPI", "GCANNNNNTGC"},
    {"BstBI", "TTCGAA"}, {"BstXI", "CCANNNNNNTGG"}, {"DpnI", "GATC"},
    {"FseI", "GGCCGGCC"}, {"HindII", "GTYRAC"}, {"MfeI", "CAATTG"},
    {"NheI", "GCTAGC"}, {"PacI", "TTAATTAA"}, {"SfiI", "GGCCNNNNNGGCC"},
    {"SnaBI", "TACGTA"}, {"XmnI", "GAANNNNTTC"}, {"ZraI", "GACGTC"}
};

// initialize
char* sequence = NULL;
size_t seq_capacity = 0;
size_t seq_length = 0;
int seq_loaded = 0;

// general base pair representation
void to_upper(char* str) {
    while (*str) {
        *str = toupper(*str);
        str++;
    }
}

// define load function
void load_fasta(const char* filename) {
    FILE* fp = fopen(filename, "r");
    if (!fp) {
        printf("[!] Failed to open %s\n", filename);
        return;
    }

    char line[MAX_LINE];
    if (sequence) {
        free(sequence); // for unlimited size (input file)
        sequence = NULL;
    }
    seq_capacity = 1024;
    seq_length = 0;
    sequence = malloc(seq_capacity);
    if (!sequence) {
        printf("[!] Memory allocation failed\n");
        fclose(fp);
        return;
    }

    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '>') continue;
        line[strcspn(line, "\r\n")] = 0;
        size_t line_len = strlen(line);

        if (seq_length + line_len + 1 > seq_capacity) {
            seq_capacity = (seq_length + line_len + 1) * 2;
            char* new_seq = realloc(sequence, seq_capacity);
            if (!new_seq) {
                printf("[!] Memory reallocation failed\n");
                free(sequence); // free malloc
                sequence = NULL;
                fclose(fp);
                return;
            }
            sequence = new_seq;
        }
        memcpy(sequence + seq_length, line, line_len);
        seq_length += line_len;
    }
    sequence[seq_length] = '\0';
    fclose(fp);
    to_upper(sequence);
    seq_loaded = 1;
    printf("[+] Loaded sequence: %lu base pairs\n", seq_length);
}

// define list function
void list_enzymes() {
    printf("Available enzymes:\n");
    for (int i = 0; i < ENZYME_COUNT; ++i) {
        printf("  %-10s %s\n", enzymes[i].name, enzymes[i].pattern);
    }
}

// define scan function
void scan_enzyme(const char* name) {
    if (!seq_loaded) {
        printf("[!] No sequence loaded. Use 'load <filename>'\n");
        return;
    }

    const char* seq = sequence;
    int found = 0;

    for (int i = 0; i < ENZYME_COUNT; ++i) {
        if (strcasecmp(name, enzymes[i].name) == 0) {
            const char* pattern = enzymes[i].pattern;
            int len = seq_length;
            int patlen = strlen(pattern);

            printf("%s (%s) found at positions:\n", enzymes[i].name, pattern);
            for (int j = 0; j <= len - patlen; ++j) {
                if (strncmp(&seq[j], pattern, patlen) == 0) {
                    printf("  [%d] ", j + 1);
                    for (int k = j - 5; k < j + patlen + 5; ++k) {
                        if (k >= 0 && k < len) {
                            if (k == j) printf("\033[1;31m"); // red start
                            putchar(seq[k]);
                            if (k == j + patlen - 1) printf("\033[0m"); // red end
                        }
                    }
                    printf("\n");
                    found++;
                }
            }
            if (!found) printf("  None found.\n");
            return;
        }
    }
    printf("[!] Enzyme '%s' not found.\n", name);
}

// help etc
void show_help() {
    printf("Commands:\n");
    printf("  load <filename>  - Load a FASTA or FSA file\n");
    printf("  list             - List available enzymes\n");
    printf("  scan <enzyme>    - Scan sequence for enzyme sites\n");
    printf("  help             - Show this help message\n");
    printf("  exit             - Exit the program\n");
}

// basic shell
int main() {
    char line[MAX_LINE];
    printf("Welcome to DNA Enzyme Shell\nType 'help' for commands.\n");

    while (1) {
        printf("Camel> ");
        if (!fgets(line, sizeof(line), stdin)) break;

        line[strcspn(line, "\r\n")] = 0;

        if (strncmp(line, "load ", 5) == 0) {
            load_fasta(line + 5);
        } else if (strcmp(line, "list") == 0) {
            list_enzymes();
        } else if (strncmp(line, "scan ", 5) == 0) {
            scan_enzyme(line + 5);
        } else if (strcmp(line, "help") == 0) {
            show_help();
        } else if (strcmp(line, "exit") == 0) {
            break;
        } else {
            printf("[!] Unknown command. Type 'help' for a list.\n");
        }
    }

    if (sequence) free(sequence);
    printf("Goodbye.\n");
    return 0;
}
