#ifndef PCHFOLD_H_
#define PCHFOLD_H_
#include <psrxml.h>
#include <phcx.h>
#ifndef uint64
#define uint64 unsigned long long int
#endif
#ifndef SPEED_OF_LIGHT
#define SPEED_OF_LIGHT 2.99792458e8
#endif

/*
 * STRUCT declarations. 
 * 
 * 
 */

typedef struct pch_tune_state {
    psrxml* header;

    double period, acc, jerk, dm;
    double period_step, period_range;
    double acc_step, acc_range;
    double jerk_step, jerk_range;
    double dm_step, dm_range;
    char use_acc, use_jerk;
    int nsubints, nchans, profile_nbins;
    int rms_size;
    int file_num;
    char* data_file_name;
} pch_tune_state_t;

typedef struct pch_tune_scrunched {
    float*** dataCube; // ->[subint,subband,bin]
    double rms;
} pch_tune_scrunched_t;

typedef struct pch_tune_stack_result {
    double best_snr;
    double best_width;
    double best_period;
    double best_accn;
    double best_jerk;
    int best_p_trial;
    int best_acc_trial;
    int best_jerk_trial;
    double best_pfactor, best_pdfactor, best_pddfactor;

    float* best_profile;
} pch_tune_stack_result_t;

typedef struct pch_tune_optimise_result {
    pch_tune_stack_result** stack_results;
    int ndms;
    int best_dm_step;
    double best_dm;
    float** best_subints;
    float** orig_subints;
    float** best_subbands;
    float** orig_subbands;
    float* orig_profile;
    phcx_snr_block* snr_block;


} pch_tune_optimise_result_t;

/*
 * FUNCTION prototypes 
 * 
 */

pch_tune_scrunched_t* pch_tune_make_scrunched(pch_tune_state_t* state);
pch_tune_optimise_result_t* pch_tune_optimise_scrunched(pch_tune_state_t* state, pch_tune_scrunched_t* scrunched);

double pch_tune_get_adjusted_time(pch_tune_state_t* state, double natural_time);
float** pch_tune_dedisperse(pch_tune_state_t* state, double dm,
        float*** profiles);
pch_tune_stack_result_t* pch_tune_stack_slide(pch_tune_state_t* state,
        float** subints, phcx_snr_block** snr_block, float rmss, int dmStep);

void pch_tune_smooth(float* pr, int maxwidth, double rmsp, int nbin,
        double* snr, double* width);
#endif /*PCHFOLD_H_*/
