#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <sys/time.h>
#include <sys/resource.h>

#define SUCCEEDED 1
#define FAILED 0
#define TRUE 1
#define FALSE 0
#define BUF_SIZE 10240
#define REF_ID_LEN_MAX 128
#define REF_SEQ_NUM_MAX 9999
#define REF_SEQ_LEN_MAX 1000000000
#define REF_SEQ_LEN_MIN 100
#define FASTQ_NUM_MAX 100000000
#define FASTQ_LEN_MAX 1000000
#define TEMPLATE_NUM_MAX 100000000
#define TEMPLATE_LEN_MAX 1000000
#define EXISTS_LINE_FEED 1
#define DEPTH_MAX 1000
#define RATIO_MAX 1000
#define PROCESS_MODEL 1
#define PROCESS_SAMPLING 2
#define PROCESS_SAMPLING_REUSE 3
#define PROCESS_SAMPLING_STORE 4
#define PROCESS_TEMPLATE 5
#define ACCURACY_MAX 100
#define STATE_MAX 50

/////////////////////////////////////////
// Definitions of structures           //
/////////////////////////////////////////

// simulation parameters
struct sim_t {
  int process;
  double depth;
  double accuracy_mean, accuracy_max, accuracy_min;
  long len_min, len_max; 
  double len_mean, len_sd; 
  long long len_quota;
  long sub_ratio, ins_ratio, del_ratio;
  double sub_rate, ins_rate, del_rate;
  int set_flg[20];
  long res_num;
  long long res_len_total; 
  double res_depth;
  double res_accuracy_mean, res_accuracy_sd;
  long res_len_min, res_len_max; 
  double res_len_mean, res_len_sd; 
  long res_sub_num, res_ins_num, res_del_num;
  double res_sub_rate, res_ins_rate, res_del_rate;
  char *prefix, *outfile_ref, *outfile_fq, *outfile_maf;
  char *hmm_model_file;
  char *profile_id, *profile_fq, *profile_stats;
  char *id_prefix;
};

// FASTQ
struct fastq_t {
  char *file;
  long num;
  long long len_total;
  long len_min, len_max;
  long num_filtered;
  long long len_total_filtered;
  long len_min_filtered, len_max_filtered;
  double len_mean_filtered, len_sd_filtered;
  double accuracy_mean_filtered, accuracy_sd_filtered;
};

// Reference
struct ref_t {
  char *file;
  char *seq;
  char id[REF_ID_LEN_MAX + 1];
  long len;
  long num_seq;
  long num;
};

// Template
struct templ_t {
  char *file;
  char *seq;
  char id[REF_ID_LEN_MAX + 1];
  long len;
  long num;
  long long len_total;
};

// Mutation
struct mut_t {
  long sub_thre[94], ins_thre[94], err_thre[94], del_thre;
  char *sub_nt_a, *sub_nt_t, *sub_nt_g, *sub_nt_c, *sub_nt_n, *ins_nt;
  char *qc, *new_qc, *tmp_qc, *seq, *new_seq, *maf_seq, *maf_ref_seq;
  long tmp_len_max;
  char seq_strand;
  long seq_left, seq_right;
};

// Quality code
struct qc_t {
  char character;
  double prob;
};

// Quality code of HMM
struct hmm_model_t {
  double ip[STATE_MAX+1];
  double ep[STATE_MAX+1][94];
  double tp[STATE_MAX+1][STATE_MAX+1];
  int exist_hmm;
};

/////////////////////////////////////////
// Global variances                    //
/////////////////////////////////////////

FILE *fp_filtered, *fp_stats, *fp_ref, *fp_fq, *fp_maf;
struct sim_t sim;
struct fastq_t fastq;
struct ref_t ref;
struct templ_t templ;
struct mut_t mut;
struct qc_t qc[94];
struct hmm_model_t hmm_model[ACCURACY_MAX+1];

long freq_len[FASTQ_LEN_MAX + 1];
long freq_accuracy[100000 + 1];

double uni_ep[ACCURACY_MAX+1][94];

/////////////////////////////////////////
// Prototypes of functions             //
/////////////////////////////////////////

int trim(char *line);
void init_sim_res();
int set_sim_param();
int get_ref_inf();
int get_ref_seq();
int get_fastq_inf();
int get_templ_inf();
int set_hmm_model();
int set_mut();
int simulate_by_sampling();
int simulate_by_model();
int simulate_by_templ();
int mutate();
int mutate_templ();
int count_digit(long num);
void print_sim_param(int seed);
void print_fastq_stats();
void print_templ_stats();
void print_simulation_stats();
void print_help();
void revcomp(char* str);
long get_time_cpu();
long get_time();

/////////////////////////////////////////
// Main                                //
/////////////////////////////////////////

int main (int argc, char** argv) {
  char *tp, *tmp_buf;
  long len;
  long num;
  long ratio;
  long i, j, k;
  long rst1, rst2;
  long t1, t2;
  double prob, rate;

  rst1 = get_time_cpu();
  t1 = get_time();

  memset(sim.set_flg, 0, sizeof(sim.set_flg));

  unsigned int seed = (unsigned int)time(NULL);

  // Variables for Option
  int opt, option_index;
  struct option long_options[] = {
    {"sample-fastq", 1, NULL, 0},
    {"depth", 1, NULL, 0},
    {"length-mean", 1, NULL, 0},
    {"length-sd", 1, NULL, 0},
    {"length-min", 1, NULL, 0},
    {"length-max", 1, NULL, 0},
    {"accuracy-mean", 1, NULL, 0},
    {"accuracy-min", 1, NULL, 0},
    {"accuracy-max", 1, NULL, 0},
    {"difference-ratio", 1, NULL, 0},
    {"hmm_model", 1, NULL, 0},
    {"prefix", 1, NULL, 0},
    {"sample-profile-id", 1, NULL, 0},
    {"seed", 1, NULL, 0},
    {"id-prefix", 1, NULL, 0},
    {"template-fasta", 1, NULL, 0},
    {0, 0, 0, 0}
  };

  // Option parsing
  option_index = 0;
  while ((opt = getopt_long(argc, argv, "", long_options, &option_index)) != -1) {
    switch (opt) {
    case 0:
      sim.set_flg[option_index] = 1;

      switch (option_index) {
      case 0:
        if ((fastq.file = (char *)malloc(strlen(optarg) + 1)) == 0) {
          fprintf(stderr, "ERROR: Cannot allocate memory.\n");
          exit(-1);
        }
        strcpy(fastq.file, optarg);
        break;

      case 1:
        sim.depth = atof(optarg);
        if (sim.depth <= 0.0) {
          fprintf(stderr, "ERROR (depth: %s): Acceptable range is more than 0.\n", optarg);
          exit(-1);
        }
        break;

      case 2:
        sim.len_mean = atof(optarg);
        if ((sim.len_mean < 1) || (sim.len_mean > FASTQ_LEN_MAX)) {
          fprintf(stderr, "ERROR (length-mean: %s): Acceptable range is 1-%d.\n", optarg, FASTQ_LEN_MAX);
          exit(-1);
        }
        break;

      case 3:
        sim.len_sd = atof(optarg);
        if ((sim.len_sd < 0) || (sim.len_sd > FASTQ_LEN_MAX)) {
          fprintf(stderr, "ERROR (length-sd: %s): Acceptable range is 0-%d.\n",
            optarg, FASTQ_LEN_MAX);
          exit(-1);
        }
        break;

      case 4:
        if (strlen(optarg) >= 8) {
          fprintf(stderr, "ERROR (length-min: %s): Acceptable range is 1-%d.\n", optarg, FASTQ_LEN_MAX);
          exit(-1);
        }
        sim.len_min = atoi(optarg);
        if ((sim.len_min < 1) || (sim.len_min > FASTQ_LEN_MAX)) {
          fprintf(stderr, "ERROR (length-min: %s): Acceptable range is 1-%d.\n", optarg, FASTQ_LEN_MAX);
          exit(-1);
        }
        break;

      case 5:
        if (strlen(optarg) >= 8) {
          fprintf(stderr, "ERROR (length-max: %s): Acceptable range is 1-%d.\n", optarg, FASTQ_LEN_MAX);
          exit(-1);
        }
        sim.len_max = atoi(optarg);
        if ((sim.len_max < 1) || (sim.len_max > FASTQ_LEN_MAX)) {
          fprintf(stderr, "ERROR (length-max: %s): Acceptable range is 1-%d.\n", optarg, FASTQ_LEN_MAX);
          exit(-1);
        }
        break;

      case 6:
        sim.accuracy_mean = atof(optarg);
        if ((sim.accuracy_mean < 0.0) || (sim.accuracy_mean > 1.0)) {
          fprintf(stderr, "ERROR (accuracy-mean: %s): Acceptable range is 0.0-1.0.\n", optarg);
          exit(-1);
        }
        break;

      case 7:
        sim.accuracy_min = atof(optarg);
        if ((sim.accuracy_min < 0.0) || (sim.accuracy_min > 1.0)) {
          fprintf(stderr, "ERROR (accuracy-min: %s): Acceptable range is 0.0-1.0.\n", optarg);
          exit(-1);
        }
        break;

      case 8:
        sim.accuracy_max = atof(optarg);
        if ((sim.accuracy_max < 0.0) || (sim.accuracy_max > 1.0)) {
          fprintf(stderr, "ERROR (accuracy-max: %s): Acceptable range is 0.0-1.0.\n", optarg);
          exit(-1);
        }
        break;

      case 9:
        if ((tmp_buf = (char *)malloc(strlen(optarg) + 1)) == 0) {
          fprintf(stderr, "ERROR: Cannot allocate memory.\n");
          exit(-1);
        }
        strcpy(tmp_buf, optarg);
        num = 0;
        tp = strtok(tmp_buf, ":");
        while (num < 3) {
          if (tp == NULL) {
            fprintf(stderr, "ERROR (difference-ratio: %s): Format is sub:ins:del.\n", optarg);
            exit(-1);
          }
          if (strlen(tp) >= 5) {
            fprintf(stderr, "ERROR (difference-ratio: %s): Acceptable range is 0-%d.\n", optarg, RATIO_MAX);
            exit(-1);
          }
          ratio = atoi(tp);
          if ((ratio < 0) || (ratio > RATIO_MAX)) {
            fprintf(stderr, "ERROR (difference-ratio: %s): Acceptable range is 0-%d.\n", optarg, RATIO_MAX);
            exit(-1);
          }
          if (num == 0) {
            sim.sub_ratio = ratio;
          } else if (num == 1) {
            sim.ins_ratio = ratio;
          } else if (num == 2) {
            sim.del_ratio = ratio;
          }
          num ++;
          tp = strtok(NULL, ":");
        }
        break;

      case 10:
        if ((sim.hmm_model_file = (char *)malloc(strlen(optarg) + 1)) == 0) {
          fprintf(stderr, "ERROR: Cannot allocate memory.\n");
          exit(-1);
        }
        strcpy(sim.hmm_model_file, optarg);
        break;

      case 11:
        if ((sim.prefix = (char *)malloc(strlen(optarg) + 1)) == 0) {
          fprintf(stderr, "ERROR: Cannot allocate memory.\n");
          exit(-1);
        }
        strcpy(sim.prefix, optarg);
        break;

      case 12:
        if ((sim.profile_id = (char *)malloc(strlen(optarg) + 1)) == 0) {
          fprintf(stderr, "ERROR: Cannot allocate memory.\n");
          exit(-1);
        }
        strcpy(sim.profile_id, optarg);
        break;

      case 13:
        seed = (unsigned int)atoi(optarg);
        break;

      case 14:
        if ((sim.id_prefix = (char *)malloc(strlen(optarg) + 1)) == 0) {
          fprintf(stderr, "ERROR: Cannot allocate memory.\n");
          exit(-1);
        }
        strcpy(sim.id_prefix, optarg);
        break;

      case 15:
        if ((templ.file = (char *)malloc(strlen(optarg) + 1)) == 0) {
          fprintf(stderr, "ERROR: Cannot allocate memory.\n");
          exit(-1);
        }
        strcpy(templ.file, optarg);
        break;

      default:
        break;
      }

    default:
      break;
    }
  }

  srand((unsigned int)seed);

  if (!(sim.set_flg[15])) {
    if (argv[optind] == (char*)'\0') {
      print_help();
      exit(-1);
    }
  }

  // Setting of simulation parameters     
  if (set_sim_param() == FAILED) {
    exit(-1);
  }
  print_sim_param(seed);

  // Quality code to error probability
  for (i=0; i<=93; i++) {
    qc[i].prob = pow(10, (double)i / -10);
    qc[i].character = (char)(i+33);
  }

  // Uniform error probability
  for (i=0; i<=ACCURACY_MAX; i++) {
    for (j=0; j<=93; j++) {
      uni_ep[i][j] = 0;
    }
    if (i == ACCURACY_MAX) {
      uni_ep[i][93] = 1.0;
      continue;
    }
    prob = 1.0 - i / 100.0;
    for (j=0; j<=93; j++) {
      if (prob == qc[j].prob) {
        uni_ep[i][j] = 1.0;
        break;
      } else if (prob > qc[j].prob) {
        rate = (prob - qc[j].prob) / (qc[j-1].prob - qc[j].prob);
        uni_ep[i][j-1] = rate;
        uni_ep[i][j] = 1 - rate;
        break;
      }
    }
  }

  // FASTQ
  if (sim.process == PROCESS_SAMPLING) {
    if ((fp_filtered = tmpfile()) == NULL) {
      fprintf(stderr, "ERROR: Cannot open temporary file\n");
      return FAILED;
    }

    if (get_fastq_inf() == FAILED) {
      exit(-1);
    }
    print_fastq_stats();
  } else if (sim.process == PROCESS_SAMPLING_STORE) {
    if ((fp_filtered = fopen(sim.profile_fq, "w+")) == NULL) {
      fprintf(stderr, "ERROR: Cannot open sample_profile\n");
      return FAILED;
    }
    if ((fp_stats = fopen(sim.profile_stats, "w+")) == NULL) {
      fprintf(stderr, "ERROR: Cannot open sample_profile\n");
      return FAILED;
    }

    if (get_fastq_inf() == FAILED) {
      exit(-1);
    }
    print_fastq_stats();
  } else if (sim.process == PROCESS_SAMPLING_REUSE) {
    if ((fp_filtered = fopen(sim.profile_fq, "r")) == NULL) {
      fprintf(stderr, "ERROR: Cannot open sample_profile\n");
      return FAILED;
    }
    if ((fp_stats = fopen(sim.profile_stats, "r")) == NULL) {
      fprintf(stderr, "ERROR: Cannot open sample_profile\n");
      return FAILED;
    }

    if (get_fastq_inf() == FAILED) {
      exit(-1);
    }
    print_fastq_stats();
  }

  // Quality code model
  if ((sim.process == PROCESS_MODEL) || (sim.process == PROCESS_TEMPLATE)) {
    if (set_hmm_model() == FAILED) {
      exit(-1);
    }
  }
  if(0){
  for(i=0; i<=ACCURACY_MAX; i++) {
    for(j=1; j<=STATE_MAX; j++) {
      for(k=0; k<=93; k++) {
        fprintf(stderr,"##### %ld : %ld : %ld : %f\n", i, j, k, hmm_model[i].ep[j][k]);
      }
    }
  }
  }
  if(0){
  for(i=0; i<=ACCURACY_MAX; i++) {
    for(j=1; j<=STATE_MAX; j++) {
      for(k=1; k<=STATE_MAX; k++) {
        fprintf(stderr,"##### %ld : %ld : %ld : %f\n", i, j, k, hmm_model[i].tp[j][k]);
      }
    }
  }
  }

  // Set mutation parameters and varianeces
  if (set_mut() == FAILED) {
    exit(-1);
  }

  // Creating simulated reads
  if (sim.process == PROCESS_TEMPLATE) {
    if (get_templ_inf() == FAILED) {
      exit(-1);
    }
    print_templ_stats();

    init_sim_res();

    sprintf(sim.outfile_fq, "%s.fastq", sim.prefix);
    if ((fp_fq = fopen(sim.outfile_fq, "w")) == NULL) {
      fprintf(stderr, "ERROR: Cannot open output file: %s\n", sim.outfile_fq);
      return FAILED;
    }

    sprintf(sim.outfile_maf, "%s.maf", sim.prefix);
    if ((fp_maf = fopen(sim.outfile_maf, "w")) == NULL) {
      fprintf(stderr, "ERROR: Cannot open output file: %s\n", sim.outfile_maf);
      return FAILED;
    }

    if (simulate_by_templ() == FAILED) {
      exit(-1);
    }

    print_simulation_stats();

    fclose(fp_fq);
    fclose(fp_maf);

  } else {
    if ((ref.file = (char *)malloc(strlen(argv[optind]) + 1)) == 0) {
      fprintf(stderr, "ERROR: Cannot allocate memory.\n");
      exit(-1);
    }
    strcpy(ref.file, argv[optind]);

    if (get_ref_inf() == FAILED) {
      exit(-1);
    }

    for (ref.num=1; ref.num<=ref.num_seq; ref.num++) {
      if (get_ref_seq() == FAILED) {
        exit(-1);
      }

      init_sim_res();

      sprintf(sim.outfile_fq, "%s_%04ld.fastq", sim.prefix, ref.num);
      if ((fp_fq = fopen(sim.outfile_fq, "w")) == NULL) {
        fprintf(stderr, "ERROR: Cannot open output file: %s\n", sim.outfile_fq);
        return FAILED;
      }

      sprintf(sim.outfile_maf, "%s_%04ld.maf", sim.prefix, ref.num);
      if ((fp_maf = fopen(sim.outfile_maf, "w")) == NULL) {
        fprintf(stderr, "ERROR: Cannot open output file: %s\n", sim.outfile_maf);
        return FAILED;
      }

      sim.len_quota = (long long)(sim.depth * ref.len);

      if (sim.process == PROCESS_MODEL) {
        if (simulate_by_model() == FAILED) {
          exit(-1);
        }
      } else {
        if (simulate_by_sampling() == FAILED) {
          exit(-1);
        }
      }

      print_simulation_stats();

      fclose(fp_fq);
      fclose(fp_maf);
    }
  }

  if ((sim.process == PROCESS_SAMPLING_STORE) || (sim.process == PROCESS_SAMPLING_REUSE)) {
    fclose(fp_filtered);
    fclose(fp_stats);
  }

  rst2 = get_time_cpu();
  t2 = get_time();

  fprintf(stderr, ":::: System utilization ::::\n\n");
  fprintf(stderr, "CPU time(s) : %ld\n", rst2 - rst1);
  fprintf(stderr, "Elapsed time(s) : %ld\n", t2 - t1);

  return(0);
}

///////////////////////////////////////
// Function: trim - Remove "\n"      //
///////////////////////////////////////

int trim(char *line) {
  int end_pos = strlen(line) - 1;

  if (line[end_pos] == '\n') {
    line[end_pos] = '\0';
    return 1;
  }
  return 0;
}

///////////////////////////////////////////////////////
// Function: get_ref_inf - Get reference information //
///////////////////////////////////////////////////////

int get_ref_inf() {
  FILE *fp;
  char line[BUF_SIZE];
  int ret;
  long max_len = 0;

  fprintf(stderr, ":::: Reference stats ::::\n\n");
  fprintf(stderr, "file name : %s\n", ref.file);
  fprintf(stderr, "\n");

  if ((fp = fopen(ref.file, "r")) == NULL) {
    fprintf(stderr, "ERROR: Cannot open file: %s\n", ref.file);
    return FAILED;
  }

  ref.num_seq = 0;
  ref.len = 0;

  while (fgets(line, BUF_SIZE, fp) != NULL) {
    ret = trim(line);

    if (line[0] == '>') {
      if (ref.num_seq != 0) {
        if (ref.len < REF_SEQ_LEN_MIN) {
          fprintf(stderr, "ERROR: Reference is too short. Acceptable length >= %d.\n", REF_SEQ_LEN_MIN);
          return FAILED;
        }
        fprintf(stderr, "ref.%ld (len:%ld) : %s\n", ref.num_seq, ref.len, ref.id);
        fclose(fp_ref);
        if (ref.len > max_len) {
          max_len = ref.len;
        }
      }

      ref.num_seq ++;
      if (ref.num_seq > REF_SEQ_NUM_MAX) {
        fprintf(stderr, "ERROR: References are too many. Max number of reference is %d.\n", REF_SEQ_NUM_MAX);
        return FAILED;
      }

      strncpy(ref.id, line + 1, REF_ID_LEN_MAX);
      ref.id[REF_ID_LEN_MAX] = '\0';

      sprintf(sim.outfile_ref, "%s_%04ld.ref", sim.prefix, ref.num_seq);
      if ((fp_ref = fopen(sim.outfile_ref, "w")) == NULL) {
        fprintf(stderr, "ERROR: Cannot open output file: %s\n", sim.outfile_ref);
        return FAILED;
      }

      ref.len = 0;

      while (ret != EXISTS_LINE_FEED) {
        if (fgets(line, BUF_SIZE, fp) == NULL) {
          break;
        }
        ret = trim(line);
      }

      fprintf(fp_ref, ">%s\n", ref.id);
    } else {
      ref.len += strlen(line);

      if (ref.len > REF_SEQ_LEN_MAX) {
        fprintf(stderr, "ERROR: Reference is too long. Acceptable length <= %d.\n", REF_SEQ_LEN_MAX);
        return FAILED;
      }

      fprintf(fp_ref, "%s\n", line);
    }
  }
  fclose(fp);

  if (ref.len < REF_SEQ_LEN_MIN) {
    fprintf(stderr, "ERROR: Reference is too short. Acceptable length >= %d.\n", REF_SEQ_LEN_MIN);
    return FAILED;
  }
  fprintf(stderr, "ref.%ld (len:%ld) : %s\n", ref.num_seq, ref.len, ref.id);
  fclose(fp_ref);
  if (ref.len > max_len) {
    max_len = ref.len;
  }

  fprintf(stderr, "\n");

  if ((ref.seq = (char *)malloc(max_len + 1)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  return SUCCEEDED;
}

////////////////////////////////////////////////////
// Function: get_ref_seq - Get reference sequence //
////////////////////////////////////////////////////

int get_ref_seq() {
  FILE *fp;
  char line[BUF_SIZE];
  long offset = 0;
  long copy_size;
  int ret;

  sprintf(sim.outfile_ref, "%s_%04ld.ref", sim.prefix, ref.num);

  if ((fp = fopen(sim.outfile_ref, "r")) == NULL) {
    fprintf(stderr, "ERROR: Cannot open file: %s\n", sim.outfile_ref);
    return FAILED;
  }

  while (fgets(line, BUF_SIZE, fp) != NULL) {
    ret = trim(line);

    if (line[0] == '>') {
      while (ret != EXISTS_LINE_FEED) {
        if (fgets(line, BUF_SIZE, fp) == NULL) {
          break;
        }
        ret = trim(line);
      }
    } else {
      copy_size = strlen(line);
      memcpy(ref.seq + offset, line, copy_size);
      offset += copy_size;
    }
  }
  fclose(fp);

  ref.seq[offset] = '\0';
  ref.len = strlen(ref.seq);

  return SUCCEEDED;
}

/////////////////////////////////////////////////////
// Function: get_fastq_inf - Get FASTQ information //
/////////////////////////////////////////////////////

int get_fastq_inf() {
  FILE *fp;
  char *tp, *item;
  char line[BUF_SIZE];
  char qc_tmp[FASTQ_LEN_MAX];
  long len;
  double prob;
  double accuracy;
  double accuracy_total = 0;
  long value;
  double variance;
  long i;
  int line_num;

  for (i=0; i<=sim.len_max; i++) {
    freq_len[i] = 0;
  }
  for (i=0; i<=100000; i++) {
    freq_accuracy[i] = 0;
  }

  fastq.num = 0;
  fastq.len_min = LONG_MAX;
  fastq.len_max = 0;
  fastq.len_total = 0;
  fastq.num_filtered = 0;
  fastq.len_min_filtered = LONG_MAX;
  fastq.len_max_filtered = 0;
  fastq.len_total_filtered = 0;

  if (sim.process == PROCESS_SAMPLING_REUSE) {
    while (fgets(line, BUF_SIZE, fp_stats) != NULL) {
      trim(line);
      tp = strtok(line, "\t");
      item = tp;
      tp = strtok(NULL, "\t");

      if (strcmp(item, "num") == 0) {
        fastq.num_filtered = atol(tp);
      } else if (strcmp(item, "len_total") == 0) {
        fastq.len_total_filtered = atol(tp);
      } else if (strcmp(item, "len_min") == 0) {
        fastq.len_min_filtered = atol(tp);
      } else if (strcmp(item, "len_max") == 0) {
        fastq.len_max_filtered = atol(tp);
      } else if (strcmp(item, "len_mean") == 0) {
        fastq.len_mean_filtered = atof(tp);
      } else if (strcmp(item, "len_sd") == 0) {
        fastq.len_sd_filtered = atof(tp);
      } else if (strcmp(item, "accuracy_mean") == 0) {
        fastq.accuracy_mean_filtered = atof(tp);
      } else if (strcmp(item, "accuracy_sd") == 0) {
        fastq.accuracy_sd_filtered = atof(tp);
      }
    }
  } else {
    if ((fp = fopen(fastq.file, "r")) == NULL) {
      fprintf(stderr, "ERROR: Cannot open file: %s\n", fastq.file);
      return FAILED;
    }

    qc_tmp[0] = '\0';
    len = 0;
    line_num = 0;

    while (fgets(line, BUF_SIZE, fp) != NULL) {
      if (trim(line) == EXISTS_LINE_FEED) {
        line_num ++;

        if (line_num == 4) {
          len += strlen(line);

          if (len > FASTQ_LEN_MAX) {
            fprintf(stderr, "ERROR: fastq is too long. Max acceptable length is %d.\n", FASTQ_LEN_MAX);
            return FAILED;
          }

          fastq.num ++;
          fastq.len_total += len;

          if (fastq.num > FASTQ_NUM_MAX) {
            fprintf(stderr, "ERROR: fastq is too many. Max acceptable number is %d.\n", FASTQ_NUM_MAX);
            return FAILED;
          }

          if (len > fastq.len_max) {
            fastq.len_max = len;
          }
          if (len < fastq.len_min) {
            fastq.len_min = len;
          }

          if ((len >= sim.len_min) && (len <= sim.len_max)) {
            strcat(qc_tmp, line);
            prob = 0.0;
            for (i=0; i<len; i++) {
              prob += qc[(int)qc_tmp[i] - 33].prob;
            }
            accuracy = 1.0 - (prob / len);

            if ((accuracy >= sim.accuracy_min) && (accuracy <= sim.accuracy_max)) {
              accuracy_total += accuracy;
              fastq.num_filtered ++;
              fastq.len_total_filtered += len;

              freq_len[len] ++;
              value = (int)(accuracy * 100000 + 0.5); 
              freq_accuracy[value] ++;

              fprintf(fp_filtered, "%s\n", qc_tmp);

              if (len > fastq.len_max_filtered) {
                fastq.len_max_filtered = len;
              }
              if (len < fastq.len_min_filtered) {
                fastq.len_min_filtered = len;
              }
            }
          }

          line_num = 0;
          qc_tmp[0] = '\0';
          len = 0;
        }
      } else {
        if (line_num == 3) {
          len += strlen(line);
          if (len > FASTQ_LEN_MAX) {
            fprintf(stderr, "ERROR: fastq is too long. Max acceptable length is %d.\n", FASTQ_LEN_MAX);
            return FAILED;
          }
          strcat(qc_tmp, line);
        }
      }
    }

    fclose(fp);

    if (fastq.num_filtered < 1) {
      fprintf(stderr, "ERROR: there is no sample-fastq in the valid range of length and accuracy.\n");
      return FAILED;
    }

    fastq.len_mean_filtered = (double)fastq.len_total_filtered / fastq.num_filtered;
    fastq.accuracy_mean_filtered = accuracy_total / fastq.num_filtered;

    variance = 0.0;
    for (i=0; i<=sim.len_max; i++) {
      if (freq_len[i] > 0) { 
        variance += pow((fastq.len_mean_filtered - i), 2) * freq_len[i];
      }
    }
    fastq.len_sd_filtered = sqrt(variance / fastq.num_filtered);

    variance = 0.0;
    for (i=0; i<=100000; i++) {
      if (freq_accuracy[i] > 0) { 
        variance += pow((fastq.accuracy_mean_filtered - i * 0.00001), 2) * freq_accuracy[i];
      }
    }
    fastq.accuracy_sd_filtered = sqrt(variance / fastq.num_filtered);

    if (sim.process == PROCESS_SAMPLING_STORE) {
      fprintf(fp_stats, "num\t%ld\n", fastq.num_filtered);
      fprintf(fp_stats, "len_total\t%lld\n", fastq.len_total_filtered);
      fprintf(fp_stats, "len_min\t%ld\n", fastq.len_min_filtered);
      fprintf(fp_stats, "len_max\t%ld\n", fastq.len_max_filtered);
      fprintf(fp_stats, "len_mean\t%f\n", fastq.len_mean_filtered);
      fprintf(fp_stats, "len_sd\t%f\n", fastq.len_sd_filtered);
      fprintf(fp_stats, "accuracy_mean\t%f\n", fastq.accuracy_mean_filtered);
      fprintf(fp_stats, "accuracy_sd\t%f\n", fastq.accuracy_sd_filtered);
    }
  }

  return SUCCEEDED;
}

/////////////////////////////////////////////////////
// Function: print_fastq_stats - Print FASTQ stats //
/////////////////////////////////////////////////////

void print_fastq_stats() {
  fprintf(stderr, ":::: FASTQ stats ::::\n\n");

  if (sim.process == PROCESS_SAMPLING_REUSE) {
    fprintf(stderr, "file name : %s\n", sim.profile_fq);
  } else {
    fprintf(stderr, "file name : %s\n", fastq.file);
    fprintf(stderr, "\n:: all reads ::\n");
    fprintf(stderr, "read num. : %ld\n", fastq.num);
    fprintf(stderr, "read total length : %lld\n", fastq.len_total);
    fprintf(stderr, "read min length : %ld\n", fastq.len_min);
    fprintf(stderr, "read max length : %ld\n", fastq.len_max);
  }

  fprintf(stderr, "\n:: filtered reads ::\n");
  fprintf(stderr, "read num. : %ld\n", fastq.num_filtered);
  fprintf(stderr, "read total length : %lld\n", fastq.len_total_filtered);
  fprintf(stderr, "read min length : %ld\n", fastq.len_min_filtered);
  fprintf(stderr, "read max length : %ld\n", fastq.len_max_filtered);
  fprintf(stderr, "read length mean (SD) : %f (%f)\n",
    fastq.len_mean_filtered, fastq.len_sd_filtered);
  fprintf(stderr, "read accuracy mean (SD) : %f (%f)\n",
    fastq.accuracy_mean_filtered, fastq.accuracy_sd_filtered);
  fprintf(stderr, "\n");
}

////////////////////////////////////////////////////////
// Function: get_templ_inf - Get Template information //
////////////////////////////////////////////////////////

int get_templ_inf() {
  FILE *fp;
  char line[BUF_SIZE];
  long seqlen, len;
  int ret;

  templ.num = 0;
  templ.len_total = 0;

  if ((fp = fopen(templ.file, "r")) == NULL) {
    fprintf(stderr, "ERROR: Cannot open file: %s\n", templ.file);
    return FAILED;
  }

  while (fgets(line, BUF_SIZE, fp) != NULL) {
    ret = trim(line);

    if (line[0] == '>') {
      templ.num ++;
      if (templ.num > TEMPLATE_NUM_MAX) {
        fprintf(stderr, "ERROR: template is too many. Max acceptable number is %d.\n", TEMPLATE_NUM_MAX);
        return FAILED;
      }
      seqlen = 0;
      while (ret != EXISTS_LINE_FEED) {
        if (fgets(line, BUF_SIZE, fp) == NULL) {
          break;
        }
        ret = trim(line);
      }
    } else {
      len = strlen(line);
      seqlen += len;
      templ.len_total += len;
      if (seqlen > TEMPLATE_LEN_MAX) {
        fprintf(stderr, "ERROR: template is too long. Max acceptable length is %d.\n", TEMPLATE_LEN_MAX);
        return FAILED;
      }
    }
  }

  if ((templ.seq = (char *)malloc(TEMPLATE_LEN_MAX + 1)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  return SUCCEEDED;
}

////////////////////////////////////////////////////////
// Function: print_templ_stats - Print Template stats //
////////////////////////////////////////////////////////

void print_templ_stats() {
  fprintf(stderr, ":::: Template stats ::::\n\n");

  fprintf(stderr, "file name : %s\n", templ.file);
  fprintf(stderr, "template num. : %ld\n", templ.num);
  fprintf(stderr, "template total length : %lld\n", templ.len_total);
  fprintf(stderr, "\n");
}

//////////////////////////////////////////////////////////
// Function: init_sim_res - Initiate simulation results //
//////////////////////////////////////////////////////////

void init_sim_res() {
  sim.res_num = 0;
  sim.res_len_total = 0;
  sim.res_sub_num = 0;
  sim.res_ins_num = 0;
  sim.res_del_num = 0;
  sim.res_len_min = LONG_MAX;
  sim.res_len_max = 0;
}

/////////////////////////////////////////////////////////
// Function: set_sim_param - Set simulation parameters //
/////////////////////////////////////////////////////////

int set_sim_param() {
  FILE *fp;
  long sum;

  // depth
  if (!(sim.set_flg[1])) {
    sim.depth = 20.0;
  }

  // length-mean
  if (!(sim.set_flg[2])) {
    sim.len_mean = 9000;
  }

  // length-sd
  if (!(sim.set_flg[3])) {
    sim.len_sd = 7000;
  }

  // length-min
  if (!(sim.set_flg[4])) {
    sim.len_min = 100;
  }

  // length-max
  if (!(sim.set_flg[5])) {
    sim.len_max = 1000000;
  }

  // accuracy-mean
  if (sim.set_flg[6]) {
    sim.accuracy_mean = int(sim.accuracy_mean * 100) * 0.01;
  } else {
    sim.accuracy_mean = 0.85;
  }

  // accuracy-min
  if (sim.set_flg[7]) {
    sim.accuracy_min = int(sim.accuracy_min * 100) * 0.01;
  } else {
    sim.accuracy_min = 0.75;
  }

  // accuracy-max
  if (sim.set_flg[8]) {
    sim.accuracy_max = int(sim.accuracy_max * 100) * 0.01;
  } else {
    sim.accuracy_max = 1.0;
  }

  // difference-ratio
  if (!(sim.set_flg[9])) {
    sim.sub_ratio = 6;
    sim.ins_ratio = 54;
    sim.del_ratio = 40;
  }

  sum = sim.sub_ratio + sim.ins_ratio + sim.del_ratio;
  sim.sub_rate = (double)sim.sub_ratio / sum;
  sim.ins_rate = (double)sim.ins_ratio / sum;
  sim.del_rate = (double)sim.del_ratio / sum;

  // prefix and outfile
  if (!(sim.set_flg[11])) {
    if ((sim.prefix = (char *)malloc(3)) == 0) {
      fprintf(stderr, "ERROR: Cannot allocate memory.\n");
      exit(-1);
    }
    strcpy(sim.prefix, "sd");
  }

  if ((sim.outfile_ref = (char *)malloc(strlen(sim.prefix) + 10)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  if ((sim.outfile_fq = (char *)malloc(strlen(sim.prefix) + 12)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  if ((sim.outfile_maf = (char *)malloc(strlen(sim.prefix) + 10)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  // prefix of read ID
  if (!(sim.set_flg[14])) {
    if ((sim.id_prefix = (char *)malloc(2)) == 0) {
      fprintf(stderr, "ERROR: Cannot allocate memory.\n");
      exit(-1);
    }
    strcpy(sim.id_prefix, "S");
  }

  // profile
  if (sim.set_flg[12]) {
    if ((sim.profile_fq = (char *)malloc(strlen(sim.profile_id) + 22)) == 0) {
      fprintf(stderr, "ERROR: Cannot allocate memory.\n");
      return FAILED;
    }

    if ((sim.profile_stats = (char *)malloc(strlen(sim.profile_id) + 22)) == 0) {
      fprintf(stderr, "ERROR: Cannot allocate memory.\n");
      return FAILED;
    }

    sprintf(sim.profile_fq, "sample_profile_%s.fastq", sim.profile_id);
    sprintf(sim.profile_stats, "sample_profile_%s.stats", sim.profile_id);
  }

  // length and accuracy
  if (sim.len_min > sim.len_max) {
    fprintf(stderr, "ERROR: length min(%ld) is greater than max(%ld).\n", sim.len_min, sim.len_max);
    return FAILED;
  }
  //if (sim.accuracy_min > sim.accuracy_mean) {
  //  fprintf(stderr, "ERROR: accuracy min(%f) is greater than mean(%f).\n", sim.accuracy_min, sim.accuracy_mean);
  //  return FAILED;
  //}

  // process
  if (sim.set_flg[15]) {
    if (!(sim.set_flg[10])) {
      fprintf(stderr, "ERROR: when using --template-fasta, --hmm_model(model-based) must be set.\n");
      return FAILED;
    }
    sim.process = PROCESS_TEMPLATE;
  } else if (sim.set_flg[10]) {
    if ((sim.set_flg[0]) || (sim.set_flg[12])) {
      fprintf(stderr, "ERROR: either --sample-fastq(and/or --sample-profile-id)(sampling-based) or --hmm_model(model-based) must be set.\n");
      return FAILED;
    }
    sim.process = PROCESS_MODEL;
  } else {
    if (sim.set_flg[0]) {
      if (sim.set_flg[12]) {
        sim.process = PROCESS_SAMPLING_STORE;
      } else {
        sim.process = PROCESS_SAMPLING;
      }
    } else {
      if (sim.set_flg[12]) {
        sim.process = PROCESS_SAMPLING_REUSE;
      } else {
        fprintf(stderr, "ERROR: either --sample-fastq(and/or --sample-profile-id)(sampling-based) or --hmm_model(model-based) must be set.\n");
        return FAILED;
      }
    }
  }

  // sample profile
  if (sim.process == PROCESS_SAMPLING_STORE) {
    if ((fp = fopen(sim.profile_fq, "r")) != NULL) {
      fprintf(stderr, "ERROR: %s exists.\n", sim.profile_fq);
      fclose(fp);
      return FAILED;
    }
    if ((fp = fopen(sim.profile_stats, "r")) != NULL) {
      fprintf(stderr, "ERROR: %s exists.\n", sim.profile_stats);
      fclose(fp);
      return FAILED;
    }
  }

  if (sim.process == PROCESS_SAMPLING_REUSE) {
    if ((fp = fopen(sim.profile_fq, "r")) == NULL) {
      fprintf(stderr, "ERROR: %s does not exist.\n", sim.profile_fq);
      return FAILED;
    }
    fclose(fp);
    if ((fp = fopen(sim.profile_stats, "r")) == NULL) {
      fprintf(stderr, "ERROR: %s does not exist.\n", sim.profile_stats);
      return FAILED;
    }
    fclose(fp);
  }

  return SUCCEEDED;
}

////////////////////////////////////////////////////////
// Function: simulate_by_sampling - Simulate by model //
////////////////////////////////////////////////////////

int simulate_by_sampling() {
  long len;
  long long len_total = 0;
  long sampling_num, sampling_interval, sampling_value, sampling_residue;
  long num;
  long i, j;
  long index;
  long value;
  double accuracy, accuracy_total = 0.0;
  double prob, variance;
  char id[128];
  int digit_num1[4], digit_num2[4], digit_num[4];

  for (i=0; i<=sim.len_max; i++) {
    freq_len[i] = 0;
  }
  for (i=0; i<=100000; i++) {
    freq_accuracy[i] = 0;
  }

  for (i=0; i<=93; i++) {
    mut.err_thre[i] = int(qc[i].prob * 1000000 + 0.5);
    mut.sub_thre[i] = int((qc[i].prob * sim.sub_rate) * 1000000 + 0.5);
    mut.ins_thre[i] = int((qc[i].prob * (sim.sub_rate + sim.ins_rate)) * 1000000 + 0.5);
  }

  sampling_num = (long)(sim.len_quota / fastq.len_total_filtered);
  sampling_residue = sim.len_quota % fastq.len_total_filtered;
  if (sampling_residue == 0) {
    sampling_interval = 1;
  } else {
    sampling_interval = (long)((double)(fastq.len_total_filtered / sampling_residue) * 2 + 0.5);
    if (sampling_interval > (long)(fastq.num_filtered * 0.5)) {
      sampling_interval = (long)(fastq.num_filtered * 0.5);
    }
  }

  // Make simulation data
  while (len_total < sim.len_quota) {
    rewind(fp_filtered);

    sampling_value = rand() % fastq.num_filtered;
    while (fgets(mut.qc, fastq.len_max_filtered + 2, fp_filtered) != NULL) {
      if (len_total >= sim.len_quota) {
        break;
      }

      trim(mut.qc);

      if (sampling_value % sampling_interval == 0) {
        num = sampling_num + 1;
      } else {
        num = sampling_num;
      }
      sampling_value ++;

      for (i=0; i<num; i++) {
        if (len_total >= sim.len_quota) {
          break;
        }

        mut.tmp_len_max = sim.len_quota - len_total;
        if (mut.tmp_len_max < sim.len_min) {
          mut.tmp_len_max = sim.len_min;
        }

        len = strlen(mut.qc);
        prob = 0.0;
        for (j=0; j<len; j++) {
          prob += qc[(int)mut.qc[j] - 33].prob;
        }
        accuracy = 1.0 - (prob / len);
        mut.del_thre = int(((1 - accuracy) * sim.del_rate) * 1000000 + 0.5);

        if (mutate() == FAILED) {
          return FAILED;
        }

        sim.res_num ++;
        len = strlen(mut.new_seq);
        sim.res_len_total += len;
        len_total += len;
        freq_len[len] ++;

        if (len > sim.res_len_max) {
          sim.res_len_max = len;
        }
        if (len < sim.res_len_min) {
          sim.res_len_min = len;
        }

        prob = 0.0;
        for (j=0; j<len; j++) {
          prob += qc[(int)mut.new_qc[j] - 33].prob;
        }
        accuracy = 1.0 - (prob / len);
        accuracy_total += accuracy;
        value = (int)(accuracy * 100000 + 0.5);
        freq_accuracy[value] ++;

        sprintf(id, "%s%ld_%ld", sim.id_prefix, ref.num, sim.res_num);
        fprintf(fp_fq, "@%s\n%s\n+%s\n%s\n", id, mut.new_seq, id, mut.new_qc);

        digit_num1[0] = 3;
        digit_num2[0] = 1 + count_digit(sim.res_num);
        digit_num[0] = (digit_num1[0] >= digit_num2[0]) ? digit_num1[0] : digit_num2[0];

        digit_num1[1] = count_digit((mut.seq_left - 1));
        digit_num2[1] = 1;
        digit_num[1] = (digit_num1[1] >= digit_num2[1]) ? digit_num1[1] : digit_num2[1];

        digit_num1[2] = count_digit((mut.seq_right - mut.seq_left + 1));
        digit_num2[2] = count_digit(len);
        digit_num[2] = (digit_num1[2] >= digit_num2[2]) ? digit_num1[2] : digit_num2[2];

        digit_num1[3] = count_digit(ref.len);
        digit_num2[3] = count_digit(len);
        digit_num[3] = (digit_num1[3] >= digit_num2[3]) ? digit_num1[3] : digit_num2[3];

        fprintf(fp_maf, "a\ns ref"); 
        while (digit_num1[0] ++ < digit_num[0]) {
          fprintf(fp_maf, " ");
        }
        while (digit_num1[1] ++ < digit_num[1]) {
          fprintf(fp_maf, " ");
        }
        fprintf(fp_maf, " %ld", mut.seq_left - 1);
        while (digit_num1[2] ++ < digit_num[2]) {
          fprintf(fp_maf, " ");
        }
        fprintf(fp_maf, " %ld +", mut.seq_right - mut.seq_left + 1);
        while (digit_num1[3] ++ < digit_num[3]) {
          fprintf(fp_maf, " ");
        }
        fprintf(fp_maf, " %ld %s\n", ref.len, mut.maf_ref_seq);
        fprintf(fp_maf, "s %s", id); 
        while (digit_num2[0] ++ < digit_num[0]) {
          fprintf(fp_maf, " ");
        }
        while (digit_num2[1] ++ < digit_num[1]) {
          fprintf(fp_maf, " ");
        }
        fprintf(fp_maf, " %d", 0);
        while (digit_num2[2] ++ < digit_num[2]) {
          fprintf(fp_maf, " ");
        }
        fprintf(fp_maf, " %ld %c", len, mut.seq_strand);
        while (digit_num2[3] ++ < digit_num[3]) {
          fprintf(fp_maf, " ");
        }
        fprintf(fp_maf, " %ld %s\n\n", len, mut.maf_seq);
      }
    }

    sampling_num = 0;
  }

  sim.res_len_mean = (double)sim.res_len_total / sim.res_num;
  sim.res_accuracy_mean = accuracy_total / sim.res_num;

  if (sim.res_num == 1) {
    sim.res_len_sd = 0.0;
    sim.res_accuracy_sd = 0.0;
  } else {
    variance = 0.0;
    for (i=0; i<=sim.len_max; i++) {
      if (freq_len[i] > 0) {
        variance += pow((sim.res_len_mean - i), 2) * freq_len[i];
      }
    }
    sim.res_len_sd = sqrt(variance / sim.res_num);

    variance = 0.0;
    for (i=0; i<=100000; i++) {
      if (freq_accuracy[i] > 0) {
        variance += pow((sim.res_accuracy_mean - i * 0.00001), 2) * freq_accuracy[i];
      }
    }
    sim.res_accuracy_sd = sqrt(variance / sim.res_num);
  }

  return SUCCEEDED;
}

/////////////////////////////////////////////////////
// Function: simulate_by_model - Simulate by Model //
/////////////////////////////////////////////////////

int simulate_by_model() {
  long len;
  long long len_total = 0;
  long num;
  long i, j, k, l;
  double prob, mean, variance, sd;
  double kappa, theta, gamma;
  double len_prob_total, freq_total, accuracy_prob_total, qc_prob_total, value, sum;
  double accuracy_total = 0.0;
  int accuracy;
  static long prob2len[100001], prob2accuracy[100001];
  static long freq2qc[ACCURACY_MAX+1][1001];
  long state;
  static long init2state[ACCURACY_MAX+1][1001];
  static long emis2qc[ACCURACY_MAX+1][STATE_MAX+1][101];
  static long tran2state[ACCURACY_MAX+1][STATE_MAX+1][101];
  long len_rand_value, accuracy_rand_value;
  long qc_rand_value_freq[ACCURACY_MAX+1];
  long qc_rand_value_init[ACCURACY_MAX+1];
  long qc_rand_value_emis[ACCURACY_MAX+1][STATE_MAX+1];
  long qc_rand_value_tran[ACCURACY_MAX+1][STATE_MAX+1];
  long start_wk, end_wk;
  long index, pre_index;
  long accuracy_min, accuracy_max;
  char id[128];
  int digit_num1[4], digit_num2[4], digit_num[4];

  for (i=0; i<=sim.len_max; i++) {
    freq_len[i] = 0;
  }
  for (i=0; i<=100000; i++) {
    freq_accuracy[i] = 0;
  }

  for (i=0; i<=93; i++) {
    mut.err_thre[i] = int(qc[i].prob * 1000000 + 0.5);
    mut.sub_thre[i] = int((qc[i].prob * sim.sub_rate) * 1000000 + 0.5);
    mut.ins_thre[i] = int((qc[i].prob * (sim.sub_rate + sim.ins_rate)) * 1000000 + 0.5);
  }

  // length distribution
  variance = pow(sim.len_sd, 2);
  kappa = pow(sim.len_mean, 2) / variance; 
  theta = variance / sim.len_mean;
  gamma = tgamma(kappa);

  if (sim.len_sd == 0.0) {
    prob2len[1] = int(sim.len_mean + 0.5);
    len_rand_value = 1;
  } else {
    start_wk = 1; 
    len_prob_total = 0.0;
    for (i=sim.len_min; i<=sim.len_max; i++) {
      len_prob_total += pow(i, kappa-1) * exp(-1 * i / theta) / pow(theta, kappa) / gamma;
      end_wk = int(len_prob_total * 100000 + 0.5);
      if (end_wk > 100000) {
        end_wk = 100000;
      }

      for (j=start_wk; j<=end_wk; j++) {
        prob2len[j] = i;
      }

      if (end_wk >= 100000) {
        break;
      }
      start_wk = end_wk + 1;
    }
    len_rand_value = end_wk;
  }

  if (len_rand_value < 1) {
    fprintf(stderr, "ERROR: length parameters are not appropriate.\n");
    return FAILED;
  }

  // accuracy distribution
  mean = sim.accuracy_mean * 100;
  accuracy_max = floor(mean * 1.05);
  accuracy_min = floor(mean * 0.75);
  if (accuracy_max > 100) {
    accuracy_max = 100;
  }

  freq_total = 0.0;
  for (i=accuracy_min; i<=accuracy_max; i++) {
    freq_total += exp(0.22 * i);
  }
  start_wk = 1; 
  accuracy_prob_total = 0.0;
  for (i=accuracy_min; i<=accuracy_max; i++) {
    accuracy_prob_total += exp(0.22 * i) / freq_total;
    end_wk = int(accuracy_prob_total * 100000 + 0.5);
    if (end_wk > 100000) {
      end_wk = 100000;
    }

    for (j=start_wk; j<=end_wk; j++) {
      prob2accuracy[j] = i;
    }

    if (end_wk >= 100000) {
      break;
    }
    start_wk = end_wk + 1;
  }
  accuracy_rand_value = end_wk;

  if (accuracy_rand_value < 1) {
    fprintf(stderr, "ERROR: accuracy parameters are not appropriate.\n");
    return FAILED;
  }

  // quality code distribution
  for (i=accuracy_min; i<=accuracy_max; i++) {

    if (hmm_model[i].exist_hmm == 1) {
      start_wk = 1; 
      qc_prob_total = 0.0;

      for (j=1; j<=STATE_MAX; j++) {
        if (hmm_model[i].ip[j] == 0) {
          continue;
        }
        qc_prob_total += hmm_model[i].ip[j];
        end_wk = int(qc_prob_total * 100 + 0.5);
        if (end_wk > 100) {
          end_wk = 100;
        }

        for (k=start_wk; k<=end_wk; k++) {
          init2state[i][k] = j;
        }

        if (end_wk >= 100) {
          break;
        }
        start_wk = end_wk + 1;
      }
      qc_rand_value_init[i] = end_wk;

      for (j=1; j<=STATE_MAX; j++) {
        start_wk = 1; 
        qc_prob_total = 0.0;

        for (k=0; k<=93; k++) {
          if (hmm_model[i].ep[j][k] == 0) {
            continue;
          }
          qc_prob_total += hmm_model[i].ep[j][k];
          end_wk = int(qc_prob_total * 100 + 0.5);
          if (end_wk > 100) {
            end_wk = 100;
          }

          for (l=start_wk; l<=end_wk; l++) {
            emis2qc[i][j][l] = k;
          }

          if (end_wk >= 100) {
            break;
          }
          start_wk = end_wk + 1;
        }
        qc_rand_value_emis[i][j] = end_wk;
      }

      for (j=1; j<=STATE_MAX; j++) {
        start_wk = 1; 
        qc_prob_total = 0.0;

        for (k=1; k<=STATE_MAX; k++) {
          if (hmm_model[i].tp[j][k] == 0) {
            continue;
          }
          qc_prob_total += hmm_model[i].tp[j][k];
          end_wk = int(qc_prob_total * 100 + 0.5);
          if (end_wk > 100) {
            end_wk = 100;
          }

          for (l=start_wk; l<=end_wk; l++) {
            tran2state[i][j][l] = k;
          }

          if (end_wk >= 100) {
            break;
          }
          start_wk = end_wk + 1;
        }
        qc_rand_value_tran[i][j] = end_wk;
      }
    } else {
      start_wk = 1; 
      qc_prob_total = 0.0;

      for (j=0; j<=93; j++) {
        if (uni_ep[i][j] == 0) {
          continue;
        }
        qc_prob_total += uni_ep[i][j];
        end_wk = int(qc_prob_total * 1000 + 0.5);
        if (end_wk > 1000) {
          end_wk = 1000;
        }

        for (k=start_wk; k<=end_wk; k++) {
          freq2qc[i][k] = j;
        }

        if (end_wk >= 1000) {
          break;
        }
        start_wk = end_wk + 1;
      }
      qc_rand_value_freq[i] = end_wk;
    }
  }

  // simulation
  while (len_total < sim.len_quota) {
    index = rand() % len_rand_value + 1;
    len = prob2len[index];
    if (len_total + len > sim.len_quota) {
      len = sim.len_quota - len_total;

      if (len < sim.len_min) {
        len = sim.len_min;
      }
    }

    mut.tmp_len_max = len;

    index = rand() % accuracy_rand_value + 1;
    accuracy = prob2accuracy[index];

    if (hmm_model[accuracy].exist_hmm == 1) {
      num = 0;
      index = rand() % qc_rand_value_init[accuracy] + 1;
      state = init2state[accuracy][index];
      index = rand() % qc_rand_value_emis[accuracy][state] + 1;
      index = emis2qc[accuracy][state][index];
      mut.qc[num ++] = qc[index].character;
      while (num < len) {
        index = rand() % qc_rand_value_tran[accuracy][state] + 1;
        state = tran2state[accuracy][state][index];
        index = rand() % qc_rand_value_emis[accuracy][state] + 1;
        index = emis2qc[accuracy][state][index];
        mut.qc[num ++] = qc[index].character;
      }
      mut.qc[num] = '\0';
    } else {
      num = 0;
      while (num < len) {
        index = rand() % qc_rand_value_freq[accuracy] + 1;
        index = freq2qc[accuracy][index];
        mut.qc[num ++] = qc[index].character;
      }
      mut.qc[num] = '\0';
    }

    mut.del_thre = int(((100 - accuracy) * sim.del_rate) * 10000 + 0.5);
    if (mutate() == FAILED) {
      return FAILED;
    }

    len = strlen(mut.new_seq);
    sim.res_len_total += len;
    len_total += len;
    freq_len[len] ++;
    sim.res_num ++;

    if (len > sim.res_len_max) {
      sim.res_len_max = len;
    }
    if (len < sim.res_len_min) {
      sim.res_len_min = len;
    }

    prob = 0.0;
    for (i=0; i<len; i++) {
      prob += qc[(int)mut.new_qc[i] - 33].prob;
    }
    value = 1.0 - (prob / len);
    accuracy_total += value;
    accuracy = (int)(value * 100000 + 0.5);
    freq_accuracy[accuracy] ++;

    sprintf(id, "%s%ld_%ld", sim.id_prefix, ref.num, sim.res_num);
    fprintf(fp_fq, "@%s\n%s\n+%s\n%s\n", id, mut.new_seq, id, mut.new_qc);

    digit_num1[0] = 3;
    digit_num2[0] = 1 + count_digit(sim.res_num);
    digit_num[0] = (digit_num1[0] >= digit_num2[0]) ? digit_num1[0] : digit_num2[0];

    digit_num1[1] = count_digit((mut.seq_left - 1));
    digit_num2[1] = 1;
    digit_num[1] = (digit_num1[1] >= digit_num2[1]) ? digit_num1[1] : digit_num2[1];

    digit_num1[2] = count_digit((mut.seq_right - mut.seq_left + 1));
    digit_num2[2] = count_digit(len);
    digit_num[2] = (digit_num1[2] >= digit_num2[2]) ? digit_num1[2] : digit_num2[2];

    digit_num1[3] = count_digit(ref.len);
    digit_num2[3] = count_digit(len);
    digit_num[3] = (digit_num1[3] >= digit_num2[3]) ? digit_num1[3] : digit_num2[3];

    fprintf(fp_maf, "a\ns ref");
    while (digit_num1[0] ++ < digit_num[0]) {
      fprintf(fp_maf, " ");
    }
    while (digit_num1[1] ++ < digit_num[1]) {
      fprintf(fp_maf, " ");
    }
    fprintf(fp_maf, " %ld", mut.seq_left - 1);
    while (digit_num1[2] ++ < digit_num[2]) {
      fprintf(fp_maf, " ");
    }
    fprintf(fp_maf, " %ld +", mut.seq_right - mut.seq_left + 1);
    while (digit_num1[3] ++ < digit_num[3]) {
      fprintf(fp_maf, " ");
    }
    fprintf(fp_maf, " %ld %s\n", ref.len, mut.maf_ref_seq);
    fprintf(fp_maf, "s %s", id);
    while (digit_num2[0] ++ < digit_num[0]) {
      fprintf(fp_maf, " ");
    }
    while (digit_num2[1] ++ < digit_num[1]) {
      fprintf(fp_maf, " ");
    }
    fprintf(fp_maf, " %d", 0);
    while (digit_num2[2] ++ < digit_num[2]) {
      fprintf(fp_maf, " ");
    }
    fprintf(fp_maf, " %ld %c", len, mut.seq_strand);
    while (digit_num2[3] ++ < digit_num[3]) {
      fprintf(fp_maf, " ");
    }
    fprintf(fp_maf, " %ld %s\n\n", len, mut.maf_seq);
  }

  sim.res_len_mean = (double)sim.res_len_total / sim.res_num;
  sim.res_accuracy_mean = accuracy_total / sim.res_num;

  if (sim.res_num == 1) {
    sim.res_len_sd = 0.0;
    sim.res_accuracy_sd = 0.0;
  } else {
    variance = 0.0;
    for (i=0; i<=sim.len_max; i++) {
      if (freq_len[i] > 0) {
        variance += pow((sim.res_len_mean - i), 2) * freq_len[i];
      }
    }
    sim.res_len_sd = sqrt(variance / sim.res_num);

    variance = 0.0;
    for (i=0; i<=100000; i++) {
      if (freq_accuracy[i] > 0) {
        variance += pow((sim.res_accuracy_mean - i * 0.00001), 2) * freq_accuracy[i];
      }
    }
    sim.res_accuracy_sd = sqrt(variance / sim.res_num);
  }

  return SUCCEEDED;
}

////////////////////////////////////////////////////////
// Function: simulate_by_templ - Simulate by Template //
////////////////////////////////////////////////////////

int simulate_by_templ() {
  FILE *fp;
  char line[BUF_SIZE];
  long offset = 0;
  long copy_size;
  int ret;
  char *ret_pointer;
  long len;
  long long len_total = 0;
  long num;
  long i, j, k, l;
  double prob, mean, variance, sd;
  double freq_total, accuracy_prob_total, qc_prob_total, value, sum;
  double accuracy_total = 0.0;
  int accuracy;
  static long prob2accuracy[100001];
  static long freq2qc[ACCURACY_MAX+1][1001];
  long state;
  static long init2state[ACCURACY_MAX+1][1001];
  static long emis2qc[ACCURACY_MAX+1][STATE_MAX+1][101];
  static long tran2state[ACCURACY_MAX+1][STATE_MAX+1][101];
  long accuracy_rand_value;
  long qc_rand_value_freq[ACCURACY_MAX+1];
  long qc_rand_value_init[ACCURACY_MAX+1];
  long qc_rand_value_emis[ACCURACY_MAX+1][STATE_MAX+1];
  long qc_rand_value_tran[ACCURACY_MAX+1][STATE_MAX+1];
  long start_wk, end_wk;
  long index, pre_index;
  long accuracy_min, accuracy_max;
  char id[128];
  int digit_num1[4], digit_num2[4], digit_num[4];

  for (i=0; i<=100000; i++) {
    freq_accuracy[i] = 0;
  }

  for (i=0; i<=93; i++) {
    mut.err_thre[i] = int(qc[i].prob * 1000000 + 0.5);
    mut.sub_thre[i] = int((qc[i].prob * sim.sub_rate) * 1000000 + 0.5);
    mut.ins_thre[i] = int((qc[i].prob * (sim.sub_rate + sim.ins_rate)) * 1000000 + 0.5);
  }

  // accuracy distribution
  mean = sim.accuracy_mean * 100;
  accuracy_max = floor(mean * 1.05);
  accuracy_min = floor(mean * 0.75);
  if (accuracy_max > 100) {
    accuracy_max = 100;
  }

  freq_total = 0.0;
  for (i=accuracy_min; i<=accuracy_max; i++) {
    freq_total += exp(0.22 * i);
  }
  start_wk = 1; 
  accuracy_prob_total = 0.0;
  for (i=accuracy_min; i<=accuracy_max; i++) {
    accuracy_prob_total += exp(0.22 * i) / freq_total;
    end_wk = int(accuracy_prob_total * 100000 + 0.5);
    if (end_wk > 100000) {
      end_wk = 100000;
    }

    for (j=start_wk; j<=end_wk; j++) {
      prob2accuracy[j] = i;
    }

    if (end_wk >= 100000) {
      break;
    }
    start_wk = end_wk + 1;
  }
  accuracy_rand_value = end_wk;

  if (accuracy_rand_value < 1) {
    fprintf(stderr, "ERROR: accuracy parameters are not appropriate.\n");
    return FAILED;
  }

  // quality code distribution
  for (i=accuracy_min; i<=accuracy_max; i++) {

    if (hmm_model[i].exist_hmm == 1) {
      start_wk = 1; 
      qc_prob_total = 0.0;

      for (j=1; j<=STATE_MAX; j++) {
        if (hmm_model[i].ip[j] == 0) {
          continue;
        }
        qc_prob_total += hmm_model[i].ip[j];
        end_wk = int(qc_prob_total * 100 + 0.5);
        if (end_wk > 100) {
          end_wk = 100;
        }

        for (k=start_wk; k<=end_wk; k++) {
          init2state[i][k] = j;
        }

        if (end_wk >= 100) {
          break;
        }
        start_wk = end_wk + 1;
      }
      qc_rand_value_init[i] = end_wk;

      for (j=1; j<=STATE_MAX; j++) {
        start_wk = 1; 
        qc_prob_total = 0.0;

        for (k=0; k<=93; k++) {
          if (hmm_model[i].ep[j][k] == 0) {
            continue;
          }
          qc_prob_total += hmm_model[i].ep[j][k];
          end_wk = int(qc_prob_total * 100 + 0.5);
          if (end_wk > 100) {
            end_wk = 100;
          }

          for (l=start_wk; l<=end_wk; l++) {
            emis2qc[i][j][l] = k;
          }

          if (end_wk >= 100) {
            break;
          }
          start_wk = end_wk + 1;
        }
        qc_rand_value_emis[i][j] = end_wk;
      }

      for (j=1; j<=STATE_MAX; j++) {
        start_wk = 1; 
        qc_prob_total = 0.0;

        for (k=1; k<=STATE_MAX; k++) {
          if (hmm_model[i].tp[j][k] == 0) {
            continue;
          }
          qc_prob_total += hmm_model[i].tp[j][k];
          end_wk = int(qc_prob_total * 100 + 0.5);
          if (end_wk > 100) {
            end_wk = 100;
          }

          for (l=start_wk; l<=end_wk; l++) {
            tran2state[i][j][l] = k;
          }

          if (end_wk >= 100) {
            break;
          }
          start_wk = end_wk + 1;
        }
        qc_rand_value_tran[i][j] = end_wk;
      }
    } else {
      start_wk = 1; 
      qc_prob_total = 0.0;

      for (j=0; j<=93; j++) {
        if (uni_ep[i][j] == 0) {
          continue;
        }
        qc_prob_total += uni_ep[i][j];
        end_wk = int(qc_prob_total * 1000 + 0.5);
        if (end_wk > 1000) {
          end_wk = 1000;
        }

        for (k=start_wk; k<=end_wk; k++) {
          freq2qc[i][k] = j;
        }

        if (end_wk >= 1000) {
          break;
        }
        start_wk = end_wk + 1;
      }
      qc_rand_value_freq[i] = end_wk;
    }
  }

  // simulation
  if ((fp = fopen(templ.file, "r")) == NULL) {
    fprintf(stderr, "ERROR: Cannot open file: %s\n", templ.file);
    return FAILED;
  }

  offset = 0;
  while (1) {
    ret_pointer = fgets(line, BUF_SIZE, fp);
    if (((ret_pointer == NULL) || (line[0] == '>')) && (offset != 0)) {
      templ.seq[offset] = '\0';

      templ.len = strlen(templ.seq);
      len = templ.len * 2;

      index = rand() % accuracy_rand_value + 1;
      accuracy = prob2accuracy[index];

      if (hmm_model[accuracy].exist_hmm == 1) {
        num = 0;
        index = rand() % qc_rand_value_init[accuracy] + 1;
        state = init2state[accuracy][index];
        index = rand() % qc_rand_value_emis[accuracy][state] + 1;
        index = emis2qc[accuracy][state][index];
        mut.qc[num ++] = qc[index].character;
        while (num < len) {
          index = rand() % qc_rand_value_tran[accuracy][state] + 1;
          state = tran2state[accuracy][state][index];
          index = rand() % qc_rand_value_emis[accuracy][state] + 1;
          index = emis2qc[accuracy][state][index];
          mut.qc[num ++] = qc[index].character;
        }
        mut.qc[num] = '\0';
      } else {
        num = 0;
        while (num < len) {
          index = rand() % qc_rand_value_freq[accuracy] + 1;
          index = freq2qc[accuracy][index];
          mut.qc[num ++] = qc[index].character;
        }
        mut.qc[num] = '\0';
      }

      mut.del_thre = int(((100 - accuracy) * sim.del_rate) * 10000 + 0.5);
      if (mutate_templ() == FAILED) {
        return FAILED;
      }

      len = strlen(mut.new_seq);
      sim.res_len_total += len;
      len_total += len;
      freq_len[len] ++;
      sim.res_num ++;

      if (len > sim.res_len_max) {
        sim.res_len_max = len;
      }
      if (len < sim.res_len_min) {
        sim.res_len_min = len;
      }

      prob = 0.0;
      for (i=0; i<len; i++) {
        prob += qc[(int)mut.new_qc[i] - 33].prob;
      }
      value = 1.0 - (prob / len);
      accuracy_total += value;
      accuracy = (int)(value * 100000 + 0.5);
      freq_accuracy[accuracy] ++;

      sprintf(id, "%s_%ld", sim.id_prefix, sim.res_num);
      fprintf(fp_fq, "@%s\n%s\n+%s\n%s\n", id, mut.new_seq, id, mut.new_qc);

      digit_num1[0] = 3;
      digit_num2[0] = 1 + count_digit(sim.res_num);
      digit_num[0] = (digit_num1[0] >= digit_num2[0]) ? digit_num1[0] : digit_num2[0];

      digit_num1[1] = count_digit((mut.seq_left - 1));
      digit_num2[1] = 1;
      digit_num[1] = (digit_num1[1] >= digit_num2[1]) ? digit_num1[1] : digit_num2[1];

      digit_num1[2] = count_digit((mut.seq_right - mut.seq_left + 1));
      digit_num2[2] = count_digit(len);
      digit_num[2] = (digit_num1[2] >= digit_num2[2]) ? digit_num1[2] : digit_num2[2];

      digit_num1[3] = count_digit(templ.len);
      digit_num2[3] = count_digit(len);
      digit_num[3] = (digit_num1[3] >= digit_num2[3]) ? digit_num1[3] : digit_num2[3];

      fprintf(fp_maf, "a\ns %s", templ.id);
      while (digit_num1[0] ++ < digit_num[0]) {
        fprintf(fp_maf, " ");
      }
      while (digit_num1[1] ++ < digit_num[1]) {
        fprintf(fp_maf, " ");
      }
      fprintf(fp_maf, " %ld", mut.seq_left - 1);
      while (digit_num1[2] ++ < digit_num[2]) {
        fprintf(fp_maf, " ");
      }
      fprintf(fp_maf, " %ld +", mut.seq_right - mut.seq_left + 1);
      while (digit_num1[3] ++ < digit_num[3]) {
        fprintf(fp_maf, " ");
      }
      fprintf(fp_maf, " %ld %s\n", templ.len, mut.maf_ref_seq);
      fprintf(fp_maf, "s %s", id);
      while (digit_num2[0] ++ < digit_num[0]) {
        fprintf(fp_maf, " ");
      }
      while (digit_num2[1] ++ < digit_num[1]) {
        fprintf(fp_maf, " ");
      }
      fprintf(fp_maf, " %d", 0);
      while (digit_num2[2] ++ < digit_num[2]) {
        fprintf(fp_maf, " ");
      }
      fprintf(fp_maf, " %ld %c", len, mut.seq_strand);
      while (digit_num2[3] ++ < digit_num[3]) {
        fprintf(fp_maf, " ");
      }
      fprintf(fp_maf, " %ld %s\n\n", len, mut.maf_seq);
    }

    if (ret_pointer == NULL) {
      break;
    }

    ret = trim(line);

    if (line[0] == '>') {
      strncpy(templ.id, line+1, REF_ID_LEN_MAX);
      templ.id[REF_ID_LEN_MAX] = '\0';
      offset = 0;
      while (ret != EXISTS_LINE_FEED) {
        if (fgets(line, BUF_SIZE, fp) == NULL) {
          break;
        }
        ret = trim(line);
      }
    } else {
      copy_size = strlen(line);
      memcpy(templ.seq + offset, line, copy_size);
      offset += copy_size;
    }
  }
  fclose(fp);

  sim.res_len_mean = (double)sim.res_len_total / sim.res_num;
  sim.res_accuracy_mean = accuracy_total / sim.res_num;

  if (sim.res_num == 1) {
    sim.res_len_sd = 0.0;
    sim.res_accuracy_sd = 0.0;
  } else {
    variance = 0.0;
    for (i=0; i<=sim.len_max; i++) {
      if (freq_len[i] > 0) {
        variance += pow((sim.res_len_mean - i), 2) * freq_len[i];
      }
    }
    sim.res_len_sd = sqrt(variance / sim.res_num);

    variance = 0.0;
    for (i=0; i<=100000; i++) {
      if (freq_accuracy[i] > 0) {
        variance += pow((sim.res_accuracy_mean - i * 0.00001), 2) * freq_accuracy[i];
      }
    }
    sim.res_accuracy_sd = sqrt(variance / sim.res_num);
  }

  return SUCCEEDED;
}

/////////////////////////////////////////////////////////////
// Function: print_sim_param - Print simulation parameters //
/////////////////////////////////////////////////////////////

void print_sim_param(int seed) {
  fprintf(stderr, ":::: Simulation parameters :::\n\n");

  if (sim.process == PROCESS_TEMPLATE) {
    fprintf(stderr, "Error simulation for templates.\n\n");
  } else if (sim.process == PROCESS_MODEL) {
    fprintf(stderr, "Simulated by stochastic model.\n\n");
  } else {
    fprintf(stderr, "Simulated by fastq sampling.\n\n");
  }

  fprintf(stderr, "prefix : %s\n", sim.prefix);
  fprintf(stderr, "id-prefix : %s\n", sim.id_prefix);
  if (sim.set_flg[12]) {
    fprintf(stderr, "sample_profile_id : %s\n", sim.profile_id);
  }

  if (!(sim.process == PROCESS_TEMPLATE)) {
    fprintf(stderr, "depth : %lf\n", sim.depth);

    if (sim.set_flg[0]) {
      fprintf(stderr, "length-mean : (sampling FASTQ)\n");
      fprintf(stderr, "length-sd : (sampling FASTQ)\n");
    } else {
      fprintf(stderr, "length-mean : %f\n", sim.len_mean);
      fprintf(stderr, "length-sd : %f\n", sim.len_sd);
    }
    fprintf(stderr, "length-min : %ld\n", sim.len_min);
    fprintf(stderr, "length-max : %ld\n", sim.len_max);
  }

  if (sim.set_flg[0]) {
    fprintf(stderr, "accuracy-mean : (sampling FASTQ)\n");
    fprintf(stderr, "accuracy-sd : (sampling FASTQ)\n");
  } else {
    fprintf(stderr, "accuracy-mean : %f\n", sim.accuracy_mean);
  }
  fprintf(stderr, "accuracy-min : %f\n", sim.accuracy_min);
  fprintf(stderr, "accuracy-max : %f\n", sim.accuracy_max);

  fprintf(stderr, "difference-ratio : %ld:%ld:%ld\n",
    sim.sub_ratio, sim.ins_ratio, sim.del_ratio);

  fprintf(stderr, "seed : %d\n", seed);

  fprintf(stderr, "\n");
}

/////////////////////////////////////////////////////////////////
// Function: set_mut - Set mutation parameters and varianeces  //
/////////////////////////////////////////////////////////////////

int set_mut() {
  mut.sub_nt_a = (char*)"TGC";
  mut.sub_nt_t = (char*)"AGC";
  mut.sub_nt_g = (char*)"ATC";
  mut.sub_nt_c = (char*)"ATG";
  mut.sub_nt_n = (char*)"ATGC";
  mut.ins_nt = (char*)"ATGC";

  if ((mut.qc = (char *)malloc(sim.len_max + 1)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  if ((mut.new_qc = (char *)malloc(sim.len_max * 2 + 1)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  if ((mut.tmp_qc = (char *)malloc(sim.len_max * 2 + 1)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  if ((mut.seq = (char *)malloc(sim.len_max * 2 + 1)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  if ((mut.new_seq = (char *)malloc(sim.len_max * 2 + 1)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  if ((mut.maf_seq = (char *)malloc(sim.len_max * 2 + 1)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  if ((mut.maf_ref_seq = (char *)malloc(sim.len_max * 2 + 1)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  return SUCCEEDED;
}

////////////////////////////////////
// Function: mutate - Mutate read //
////////////////////////////////////

int mutate() {
  char *line;
  char nt;
  long num;
  long i, j;
  long index;
  long rand_value;
  long qc_value;
  long len;
  long offset, seq_offset, maf_offset;
  long del_num_total, pos;

  len = strlen(mut.qc);
  if (mut.tmp_len_max < len) {
    len = mut.tmp_len_max;
  }

  // Place deletions
  long del_num[len];
  for(i=0; i<len; i++) {
    del_num[i] = 0;
  }
  del_num_total = int((float)len * mut.del_thre / 1000000 + 0.5);
  num = 0;
  while (num < del_num_total) {
    pos = rand() % (len-1);
    rand_value = rand() % 1000000;
    qc_value = (int)mut.qc[pos] - 33;
    if (rand_value < mut.err_thre[qc_value]) {
      del_num[pos] ++;
      num ++;
    }
  }
  offset = 0;
  for (i=0; i<len-1; i++) {
    mut.tmp_qc[offset ++] = mut.qc[i];
    for (j=0; j<del_num[i]; j++) {
      mut.tmp_qc[offset ++] = ' ';
    }
  }
  mut.tmp_qc[offset ++] = mut.qc[len - 1];
  mut.tmp_qc[offset] = '\0';
  sim.res_del_num += del_num_total;

  len = strlen(mut.tmp_qc);

  if (len >= ref.len) {
    offset = 0;
    len = ref.len;
  } else {
    offset = rand() % (ref.len - len + 1);
  }

  mut.seq_left = offset + 1;
  mut.seq_right = offset + len;

  if (sim.res_num % 2 == 0) {
    mut.seq_strand = '+';

    for (i=0; i<len; i++) {
      nt = toupper(ref.seq[offset + i]);
      mut.seq[i] = nt;
    }
  } else {
    mut.seq_strand = '-';

    for (i=0; i<len; i++) {
      nt = toupper(ref.seq[offset + i]);

      if (nt == 'A') {
        mut.seq[len-1-i] = 'T';
      } else if (nt == 'T') {
        mut.seq[len-1-i] = 'A';
      } else if (nt == 'G') {
        mut.seq[len-1-i] = 'C';
      } else if (nt == 'C') {
        mut.seq[len-1-i] = 'G';
      } else {
        mut.seq[len-1-i] = nt;
      }
    }
  }
  mut.seq[len] = '\0';

  // Place substitutions and insertions
  offset = 0;
  seq_offset = 0;
  maf_offset = 0;
  for (i=0; i<len; i++) {
    nt = mut.seq[seq_offset ++];

    if (mut.tmp_qc[i] == ' ') {
      mut.maf_seq[maf_offset] = '-';
      mut.maf_ref_seq[maf_offset] = nt;
      maf_offset ++;
      continue;
    }

    mut.new_qc[offset] = mut.tmp_qc[i];

    rand_value = rand() % 1000000;
    qc_value = (int)mut.tmp_qc[i] - 33;

    if (rand_value < mut.sub_thre[qc_value]) {
      sim.res_sub_num ++;
      index = rand() % 3;
      if (nt == 'A') {
        mut.new_seq[offset] = mut.sub_nt_a[index];
      } else if (nt == 'T') {
        mut.new_seq[offset] = mut.sub_nt_t[index];
      } else if (nt == 'G') {
        mut.new_seq[offset] = mut.sub_nt_g[index];
      } else if (nt == 'C') {
        mut.new_seq[offset] = mut.sub_nt_c[index];
      } else {
        index = rand() % 4;
        mut.new_seq[offset] = mut.sub_nt_n[index];
      }
      mut.maf_ref_seq[maf_offset] = nt;
    } else if (rand_value < mut.ins_thre[qc_value]) {
      sim.res_ins_num ++;
      index = rand() % 8;
      if (index >= 4) {
        mut.new_seq[offset] = nt;
      } else {
        mut.new_seq[offset] = mut.ins_nt[index];
      }
      seq_offset --;
      if (mut.seq_strand == '+') {
        mut.seq_right --;
      } else {
        mut.seq_left ++;
      }
      mut.maf_ref_seq[maf_offset] = '-';
    } else {
      mut.new_seq[offset] = nt;
      mut.maf_ref_seq[maf_offset] = nt;
    }
    mut.maf_seq[maf_offset] = mut.new_seq[offset];
    maf_offset ++;
    offset ++;
  }
  mut.new_qc[offset] = '\0';
  mut.new_seq[offset] = '\0';
  mut.maf_seq[maf_offset] = '\0';
  mut.maf_ref_seq[maf_offset] = '\0';

  if (mut.seq_strand == '-') {
    revcomp(mut.maf_seq);
    revcomp(mut.maf_ref_seq);
  }

  return SUCCEEDED;
}

//////////////////////////////////////////
// Function: mutate_templ - Mutate read //
//////////////////////////////////////////

int mutate_templ() {
  char *line;
  char nt;
  long num;
  long i, j;
  long index;
  long rand_value;
  long qc_value;
  long len, seqlen;
  long offset, new_offset, seq_offset, qc_offset, maf_offset;
  long del_num_total, pos;

  len = strlen(mut.qc);

  // Place deletions
  long del_num[len];
  for(i=0; i<len; i++) {
    del_num[i] = 0;
  }
  del_num_total = int((float)len * mut.del_thre / 1000000 + 0.5);
  num = 0;
  while (num < del_num_total) {
    pos = rand() % (len-1);
    rand_value = rand() % 1000000;
    qc_value = (int)mut.qc[pos] - 33;
    if (rand_value < mut.err_thre[qc_value]) {
      del_num[pos] ++;
      num ++;
    }
  }
  offset = 0;
  for (i=0; i<len-1; i++) {
    mut.tmp_qc[offset ++] = mut.qc[i];
    for (j=0; j<del_num[i]; j++) {
      mut.tmp_qc[offset ++] = ' ';
    }
  }
  mut.tmp_qc[offset ++] = mut.qc[len - 1];
  mut.tmp_qc[offset] = '\0';

  seqlen = strlen(templ.seq);
  mut.seq_left = 1;
  mut.seq_right = seqlen;
  mut.seq_strand = '+';

  for (i=0; i<seqlen; i++) {
    nt = toupper(templ.seq[i]);
    mut.seq[i] = nt;
  }
  mut.seq[seqlen] = '\0';

  // Place substitutions and insertions
  qc_offset = 0;
  new_offset = 0;
  seq_offset = 0;
  maf_offset = 0;
  while (seq_offset < seqlen) {
    nt = mut.seq[seq_offset ++];

    if (mut.tmp_qc[qc_offset] == ' ') {
      sim.res_del_num ++;
      mut.maf_seq[maf_offset] = '-';
      mut.maf_ref_seq[maf_offset] = nt;
      maf_offset ++;
      qc_offset ++;
      continue;
    }

    mut.new_qc[new_offset] = mut.tmp_qc[qc_offset];

    rand_value = rand() % 1000000;
    qc_value = (int)mut.tmp_qc[qc_offset] - 33;

    if (rand_value < mut.sub_thre[qc_value]) {
      sim.res_sub_num ++;
      index = rand() % 3;
      if (nt == 'A') {
        mut.new_seq[new_offset] = mut.sub_nt_a[index];
      } else if (nt == 'T') {
        mut.new_seq[new_offset] = mut.sub_nt_t[index];
      } else if (nt == 'G') {
        mut.new_seq[new_offset] = mut.sub_nt_g[index];
      } else if (nt == 'C') {
        mut.new_seq[new_offset] = mut.sub_nt_c[index];
      } else {
        index = rand() % 4;
        mut.new_seq[new_offset] = mut.sub_nt_n[index];
      }
      mut.maf_ref_seq[maf_offset] = nt;
    } else if (rand_value < mut.ins_thre[qc_value]) {
      sim.res_ins_num ++;
      index = rand() % 8;
      if (index >= 4) {
        mut.new_seq[new_offset] = nt;
      } else {
        mut.new_seq[new_offset] = mut.ins_nt[index];
      }
      seq_offset --;
      mut.maf_ref_seq[maf_offset] = '-';
    } else {
      mut.new_seq[new_offset] = nt;
      mut.maf_ref_seq[maf_offset] = nt;
    }
    mut.maf_seq[maf_offset] = mut.new_seq[new_offset];
    maf_offset ++;
    qc_offset ++;
    new_offset ++;
  }
  mut.new_qc[new_offset] = '\0';
  mut.new_seq[new_offset] = '\0';
  mut.maf_seq[maf_offset] = '\0';
  mut.maf_ref_seq[maf_offset] = '\0';

  return SUCCEEDED;
}

////////////////////////////////////////////////////////////////
// Function: print_simulation_stats - Print Simulation Stats. //
////////////////////////////////////////////////////////////////

void print_simulation_stats() {
  if (sim.process == PROCESS_TEMPLATE) {
    sim.res_sub_rate = (double)sim.res_sub_num / sim.res_len_total;
    sim.res_ins_rate = (double)sim.res_ins_num / sim.res_len_total;
    sim.res_del_rate = (double)sim.res_del_num / sim.res_len_total;

    fprintf(stderr, ":::: Simulation stats (templates) ::::\n\n");
    fprintf(stderr, "read num. : %ld\n", sim.res_num);
    fprintf(stderr, "read length mean (SD) : %f (%f)\n",
      sim.res_len_mean, sim.res_len_sd);
    fprintf(stderr, "read length min : %ld\n", sim.res_len_min);
    fprintf(stderr, "read length max : %ld\n", sim.res_len_max);
    fprintf(stderr, "read accuracy mean (SD) : %f (%f)\n",
      sim.res_accuracy_mean, sim.res_accuracy_sd);
    fprintf(stderr, "substitution rate. : %f\n", sim.res_sub_rate);
    fprintf(stderr, "insertion rate. : %f\n", sim.res_ins_rate);
    fprintf(stderr, "deletion rate. : %f\n", sim.res_del_rate);
    fprintf(stderr, "\n");
  } else {
    sim.res_depth = (double)sim.res_len_total / ref.len;
    sim.res_sub_rate = (double)sim.res_sub_num / sim.res_len_total;
    sim.res_ins_rate = (double)sim.res_ins_num / sim.res_len_total;
    sim.res_del_rate = (double)sim.res_del_num / sim.res_len_total;

    fprintf(stderr, ":::: Simulation stats (ref.%ld) ::::\n\n", ref.num);
    fprintf(stderr, "read num. : %ld\n", sim.res_num);
    fprintf(stderr, "depth : %lf\n", sim.res_depth);
    fprintf(stderr, "read length mean (SD) : %f (%f)\n",
      sim.res_len_mean, sim.res_len_sd);
    fprintf(stderr, "read length min : %ld\n", sim.res_len_min);
    fprintf(stderr, "read length max : %ld\n", sim.res_len_max);
    fprintf(stderr, "read accuracy mean (SD) : %f (%f)\n",
      sim.res_accuracy_mean, sim.res_accuracy_sd);
    fprintf(stderr, "substitution rate. : %f\n", sim.res_sub_rate);
    fprintf(stderr, "insertion rate. : %f\n", sim.res_ins_rate);
    fprintf(stderr, "deletion rate. : %f\n", sim.res_del_rate);
    fprintf(stderr, "\n");
  }
}

///////////////////////////////////////////////////////
// Function: set_hmm_model - Set quality code model  //
///////////////////////////////////////////////////////

int set_hmm_model() {
  FILE *fp;
  char line[BUF_SIZE];
  char *tp;
  int accuracy, state, num;
  int i, j, k;

  if ((fp = fopen(sim.hmm_model_file, "r")) == NULL) {
    fprintf(stderr, "ERROR: Cannot open file: %s\n", sim.hmm_model_file);
    return FAILED;
  }

  for (i=0; i<=ACCURACY_MAX; i++) {
    hmm_model[i].exist_hmm = 0;
    for (j=0; j<=STATE_MAX; j++) {
      hmm_model[i].ip[j] = 0.0;
      for (k=0; k<=93; k++) {
        hmm_model[i].ep[j][k] = 0.0;
      }
      for (k=0; k<=STATE_MAX; k++) {
        hmm_model[i].tp[j][k] = 0.0;
      }
    }
  }

  while (fgets(line, BUF_SIZE, fp) != NULL) {
    trim(line);
    tp = strtok(line, " ");
    accuracy = atoi(tp);
    hmm_model[accuracy].exist_hmm = 1;

    tp = strtok(NULL, " ");
    if (strcmp(tp, "IP") == 0) {
      tp = strtok(NULL, " ");
      state = atoi(tp);
      tp = strtok(NULL, " ");
      hmm_model[accuracy].ip[state] = atof(tp);

    } else if (strcmp(tp, "EP") == 0) {
      tp = strtok(NULL, " ");
      state = atoi(tp);
      num = 0;
      tp = strtok(NULL, " ");
      while (tp != NULL) {
        hmm_model[accuracy].ep[state][num] = atof(tp);
        num ++;
        tp = strtok(NULL, " ");
      }

    } else if (strcmp(tp, "TP") == 0) {
      tp = strtok(NULL, " ");
      state = atoi(tp);
      num = 0;
      tp = strtok(NULL, " ");
      while (tp != NULL) {
        num ++;
        hmm_model[accuracy].tp[state][num] = atof(tp);
        tp = strtok(NULL, " ");
      }
    }
  }
  fclose(fp);

  return SUCCEEDED;
}

///////////////////////////////////////////////////////
// Function: get_time_cpu - Get CPU time             //
///////////////////////////////////////////////////////

long get_time_cpu() {
  struct rusage ru;
  getrusage(RUSAGE_SELF, &ru);
  return ru.ru_utime.tv_sec;
}

///////////////////////////////////////////////////////
// Function: get_time - Get time                     //
///////////////////////////////////////////////////////

long get_time() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec;
}

///////////////////////////////////////
// Function: print_help - Print help //
///////////////////////////////////////

void print_help() {
  fprintf(stderr, "\n");
  fprintf(stderr, "USAGE: pbsim [options] <reference>\n\n");
  fprintf(stderr, " <reference>           FASTA format file (text file only).\n");
  fprintf(stderr, "\n");
  fprintf(stderr, " [general options]\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  --prefix             prefix of output files (sd).\n");
  fprintf(stderr, "  --id-prefix          prefix of read ID (S).\n");
  fprintf(stderr, "  --depth              depth of coverage (20.0).\n");
  fprintf(stderr, "  --length-min         minimum length (100).\n");
  fprintf(stderr, "  --length-max         maximum length (1000000).\n");
  fprintf(stderr, "  --difference-ratio   ratio of differences (6:50:54).\n");
  fprintf(stderr, "                       (substitution:insertion:deletion)\n");
  fprintf(stderr, "                       Each value must be 0-1000, e.g. 1000:1:0 is OK.\n");
  fprintf(stderr, "                       Note that the above default value is for PacBio, and\n");
  fprintf(stderr, "                       for Nanopore we recommend increasing substitution rate\n");
  fprintf(stderr, "                       and decreasing insertion rate, such as 23:31:46.\n");
  fprintf(stderr, "  --seed               for a pseudorandom number generator (Unix time).\n");
  fprintf(stderr, "\n");
  fprintf(stderr, " [options of sampling-based simulation]\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  --sample-fastq       FASTQ format file to sample (text file only).\n");
  fprintf(stderr, "  --sample-profile-id  sample-fastq (filtered) profile ID.\n");
  fprintf(stderr, "                       when using --sample-fastq, profile is stored.\n");
  fprintf(stderr, "                       'sample_profile_<ID>.fastq', and\n");
  fprintf(stderr, "                       'sample_profile_<ID>.stats' are created.\n");
  fprintf(stderr, "                       when not using --sample-fastq, profile is re-used.\n");
  fprintf(stderr, "                       Note that when profile is used, --length-min,max,\n");
  fprintf(stderr, "                       --accuracy-min,max would be the same as the profile.\n");
  fprintf(stderr, "  --accuracy-min       minimum accuracy (0.75).\n");
  fprintf(stderr, "  --accuracy-max       maximum accuracy (1.00).\n");
  fprintf(stderr, "\n");
  fprintf(stderr, " [options of model-based simulation].\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  --hmm_model          HMM model of quality code.\n");
  fprintf(stderr, "  --length-mean        mean of length model (9000.0).\n");
  fprintf(stderr, "  --length-sd          standard deviation of length model (7000.0).\n");
  fprintf(stderr, "  --accuracy-mean      mean of accuracy model (0.85).\n");
  fprintf(stderr, "\n");
  fprintf(stderr, " [option of error simulation for templates].\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  --template-fasta     FASTA format file of templates (text file only).\n");
  fprintf(stderr, "\n");
}

/////////////////////////////////////////
// Function: count_digit - count digit //
/////////////////////////////////////////

int count_digit(long num) {
  int digit = 1;
  int quotient;

  quotient = int(num / 10);

  while (quotient != 0) {
    digit ++;
    quotient = int(quotient / 10);
  }

  return digit;
}  

//////////////////////////////////////////////////////
// Function: revcomp - convert to reverse complement//
//////////////////////////////////////////////////////

void revcomp(char* str) {
  int i, len;
  char c;

  len = strlen(str);

  for(i=0; i<len/2; i++) {
    c = str[i];
    str[i] = str[len-i-1];
    str[len-i-1] = c;
  }

  for(i=0; i<len; i++) {
    if (str[i] == 'A') {
      str[i] = 'T';
    } else if (str[i] == 'T') {
      str[i] = 'A';
    } else if (str[i] == 'G') {
      str[i] = 'C';
    } else if (str[i] == 'C') {
      str[i] = 'G';
    }
  }
}  
