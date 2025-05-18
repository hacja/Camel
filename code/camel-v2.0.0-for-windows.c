/*This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#define _CRT_SECURE_NO_WARNINGS  // close "sscanf()" security warning

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <curl/curl.h>

#define MAX_LINE 256
#define ENZYME_COUNT 55

 // Enzyme structure
typedef struct {
    char name[16];
    char pattern[32];
} Enzyme;

typedef struct {
    char* memory;
    size_t size;
} MemoryStruct;

// Enzyme data
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
    {"XbaI", "TCTAGA"}, {"XhoI", "CTCGAG"}, {"XmaI", "CCCGGG"},
    {"BsaI", "GGTCTC"}, {"BsmBI", "CGTCTC"},{"BsrGI", "TGTACA"},
    {"BstEII", "GGTNACC"}, {"EagI", "CGGCCG"}, {"AatII", "GACGTC"},
    {"BbvCI", "CCTCAGC"}, {"BsiWI", "CGTACG"},{"BspEI", "TCCGGA"},
    {"BsrI", "ACTGG"},{"BstBI", "TTCGAA"}, {"DpnI", "GATC"},
    {"FseI", "GGCCGGCC"}, {"HindII", "GTYRAC"}, {"MfeI", "CAATTG"},
    {"NheI", "GCTAGC"}, {"PacI", "TTAATTAA"}, {"SnaBI", "TACGTA"}, {"ZraI", "GACGTC"}
};

// Global variable for dna sequence (struct like)
char* sequence = NULL; // Initially NULL
size_t seq_capacity = 0; // Initial capacity
size_t seq_length = 0; // Current length
int seq_loaded = 0; // Flag for sequence loaded

// Switch strings to upper
static void to_upper(char* str) {
    while (*str) {
		*str = toupper(*str); // Convert to uppercase
		str++; // Move to next character
    }
}
// Function for download write memory
static size_t write_memory_callback(void* contents, size_t size, size_t nmemb, void* userp) {
	size_t realsize = size * nmemb; // Calculate real size
	MemoryStruct* mem = (MemoryStruct*)userp; // Cast user pointer to MemoryStruct 

	char* ptr = realloc(mem->memory, mem->size + realsize + 1); // Reallocate memory
	if (!ptr) return 0; // Check for memory allocation failure

	mem->memory = ptr; // Assign new memory to struct
	memcpy(&(mem->memory[mem->size]), contents, realsize); // Copy contents to memory
	mem->size += realsize; // Update size
	mem->memory[mem->size] = 0; // Null-terminate the string

	return realsize; // Return the size of the data
}

// Download_fasta function
static int download_fasta(const char* url, const char* out_file) {
	CURL* curl; // Initialize curl
	CURLcode res; // Initialize curl code
	FILE* fp = fopen(out_file, "wb"); // Open file for writing
    if (!fp) { 
        printf("[!] Failed to create output file.\n"); 
        return 0;
    }

	curl_global_init(CURL_GLOBAL_DEFAULT); // Initialize curl globally
	curl = curl_easy_init(); // Initialize curl
    if (!curl) {
		fclose(fp); // Close file if curl initialization fails
        return 0;
    }

	curl_easy_setopt(curl, CURLOPT_URL, url); // Set URL
	curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, NULL); // Set write function to NULL
	curl_easy_setopt(curl, CURLOPT_WRITEDATA, fp); // Set write data to file pointer

    res = curl_easy_perform(curl);
    if (res != CURLE_OK) {
		printf("[!] curl_easy_perform() failed: %s\n", curl_easy_strerror(res)); // Print error
		fclose(fp); // Close file if download fails
		curl_easy_cleanup(curl); // Cleanup curl
        return 0;
    }

    curl_easy_cleanup(curl);
    curl_global_cleanup();
	fclose(fp); // Close file
    printf("[+] Downloaded FASTA file to: %s\n", out_file);
    return 1;
}

// Load local fasta file
static void load_fasta(const char* filename) {
    FILE* fp = fopen(filename, "r");
    if (!fp) {
		printf("[!] Failed to open %s\n", filename); // Print error
        return;
    }

    char line[MAX_LINE];
    if (sequence) {
		free(sequence); // Free previous sequence
		sequence = NULL; // Initially set to NULL
    }

	seq_capacity = 1024; // Initial capacity
	seq_length = 0; // Current length
	sequence = malloc(seq_capacity); // Allocate memory for sequence
    if (!sequence) { 
		printf("[!] Memory allocation failed\n"); // Print error
        fclose(fp);
        return;
    }

    while (fgets(line, sizeof(line), fp)) {
		if (line[0] == '>') continue; // Skip header lines
		line[strcspn(line, "\r\n")] = 0; // Remove newline characters
        size_t line_len = strlen(line);

        if (seq_length + line_len + 1 > seq_capacity) {
			seq_capacity = (seq_length + line_len + 1) * 2; // Double the capacity
            char* new_seq = realloc(sequence, seq_capacity);
            if (!new_seq) {
				printf("[!] Memory reallocation failed\n"); // Print error
				free(sequence); // Free previous sequence
				sequence = NULL; // Set to NULL
				fclose(fp); // Close file
                return;
            }
			sequence = new_seq; // Assign new memory to sequence
        }

		memcpy(sequence + seq_length, line, line_len); // Copy line to sequence
		seq_length += line_len; // Update length
    }

	sequence[seq_length] = '\0'; // Null-terminate the sequence
	fclose(fp); // Close file
	to_upper(sequence); // Convert sequence to uppercase
    seq_loaded = 1;
    printf("[+] Loaded sequence: %zu base pairs\n", seq_length);
}

// Scan enzyme
static void scan_enzyme(const char* name) {
    if (!seq_loaded) {
		printf("[!] No sequence loaded. Use 'load <filename>'\n"); // Print error
        return;
    }

    for (int i = 0; i < ENZYME_COUNT; ++i) {
        if (_stricmp(name, enzymes[i].name) == 0) {
			const char* pattern = enzymes[i].pattern; // Get pattern
			int len = (int)seq_length; // Get sequence length
			int patlen = (int)strlen(pattern); // Get pattern length
			int found = 0; // Initialize found count

            printf("%s (%s) found at positions:\n", enzymes[i].name, pattern);
            for (int j = 0; j <= len - patlen; ++j) {
                if (strncmp(&sequence[j], pattern, patlen) == 0) {
                    printf("  [%d] ", j + 1);
                    for (int k = j - 5; k < j + patlen + 5; ++k) {
                        if (k >= 0 && k < len) {
                            if (k == j) printf("\033[1;31m");
                            putchar(sequence[k]);
                            if (k == j + patlen - 1) printf("\033[0m");
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

// list enzymes
static void list_enzymes() {
    printf("Available enzymes:\n");
    for (int i = 0; i < ENZYME_COUNT; ++i) {
        printf("  %-10s %s\n", enzymes[i].name, enzymes[i].pattern);
    }
}

// show help function
static void show_help() {
    printf("Commands:\n");
    printf("  download <url> <output> - Download FASTA file\n");
    printf("  load <filename>         - Load a FASTA or FSA file\n");
    printf("  list                    - List available enzymes\n");
    printf("  scan <enzyme>           - Scan sequence for enzyme sites\n");
    printf("  help                    - Show this help message\n");
    printf("  clear / cls             - Clear screen\n");
    printf("  quit / exit             - Exit the program\n");
}

// main loop
int main() {
    char line[MAX_LINE];
    printf("Welcome to DNA Enzyme Shell\nType 'help' for commands.\n");

    while (1) {
        printf("Camel> ");
        if (!fgets(line, sizeof(line), stdin)) break;
        line[strcspn(line, "\r\n")] = 0;

        if (strncmp(line, "load ", 5) == 0) {
            load_fasta(line + 5);
        }
        else if (strncmp(line, "download ", 9) == 0) {
            char url[256], out[128];
            if (sscanf(line + 9, "%255s %127s", url, out) == 2) {
                download_fasta(url, out);
            }
            else {
                printf("[!] Usage: download <url> <output_file>\n");
            }
        }
        else if (strcmp(line, "list") == 0) {
            list_enzymes();
        }
        else if (strncmp(line, "scan ", 5) == 0) {
            scan_enzyme(line + 5);
        }
        else if (strcmp(line, "help") == 0) {
            show_help();
        }
        else if (strcmp(line, "clear") == 0 || strcmp(line, "cls") == 0) {
            printf("\033[H\033[J");
        }
        else if (strcmp(line, "exit") == 0 || strcmp(line, "quit") == 0) {
            break;
        }
        else {
            printf("[!] Unknown command. Type 'help' for a list.\n");
        }
    }

    if (sequence) free(sequence);
    printf("Goodbye.\n");
    return 0;
}