#include "pch-tune.h"
#include "cbarycentre.h"
#include <phcx.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NOTSETVAL 1e200
#define NOTSETLIMIT 1e199

void pch_tune_copy_phcx_to_state(pch_tune_state_t* state, phcx* infile);

void pch_tune_print_usage() {
    printf("pch-tune [options] psrxml_file\n\n");
    printf("Options:\n\n");
    printf("--help,-h\tPrint this help message\n");

    printf("\nFolding options\n");
    printf("--period,-p [val]\tSet the folding centre period (ms,topo)\n");
    printf("--bary-period,-P [val]\tSet the folding centre period (ms,bary)\n");
    printf("--dm,-d [val]\tSet the centre dm (cm.pc^-3)\n");
    printf("--accn,-a [val]\tSet the folding acceleration (m/s/s)\n");
    printf("--jerk,-j [val]\tSet the folding jerk (m/s/s/s)\n");
    printf("--nbins,-n [val]\tSet the ouput number of profile bins\n");
    printf("--nsub,-s [val]\tSet the output numbber of sub-integrations\n");
    printf("--rms-size [val]\tSet the number of bins used to compute rms\n");


    printf("\nOptimising options\n");
    printf("--tune-accn,-A\tOptimise in acceleration\n");
    printf("--tune-jerk,-J\tOptimise in jerk\n");
    printf("--period-range [val]\tSet the period search range (ms)\n");
    printf("--period-step [val]\tSet the period search step (ms)\n");
    printf("--dm-range [val]\tSet the dm search range (cm.pc^-3)\n");
    printf("--dm-step [val]\tSet the dm search step (cm.pc^-3)\n");
    printf("--accn-range [val]\tSet the accn search range (m/s/s)\n");
    printf("--accn-step [val]\tSet the accn search step (m/s/s)\n");
    printf("--jerk-range [val]\tSet the jerk search range (m/s/s/s)\n");
    printf("--jerk-step [val]\tSet the jerk search step (m/s/s/s)\n");
    printf("\nInput options\n");
    printf("--input-phcx [file]\tReads parameters from existing phcx file\n");

    printf("\nOutput options\n");
    printf("--output-phcx [file]\tWrites pulsarhunter xml candidate files\n");
    printf("--ascii-dump\tWrites out ascii profile and subints\n");

}

int main(int argc, char** argv) {
    const char* args = "Aa:d:Jj:n:P:p:s:";
    char c;
    char ascii_dump; // writes out lots of nice ascii files for plotting
    char have_input_phcx = 0, have_output_phcx = 0, have_input_rawcand = 0,
            have_output_rawcand = 0;
    char input_phcx_filename[256], output_phcx_filename[256];
    char input_rawcand_filename[256], output_rawcand_filename[256];
    phcx* input_phcx = NULL;
    int long_opt_idx = 0;
    int opt_flag = 0;
    char dm_set = 0, convert_from_bary = 0;
    struct option long_opt[256];
    psrxml* header = (psrxml*) malloc(sizeof (psrxml));

    long_opt[long_opt_idx].name = "help";
    long_opt[long_opt_idx].has_arg = no_argument;
    long_opt[long_opt_idx].flag = NULL;
    long_opt[long_opt_idx++].val = 'h';

    long_opt[long_opt_idx].name = "period";
    long_opt[long_opt_idx].has_arg = required_argument;
    long_opt[long_opt_idx].flag = NULL;
    long_opt[long_opt_idx++].val = 'p';

    long_opt[long_opt_idx].name = "bary-period";
    long_opt[long_opt_idx].has_arg = required_argument;
    long_opt[long_opt_idx].flag = NULL;
    long_opt[long_opt_idx++].val = 'P';

    long_opt[long_opt_idx].name = "dm";
    long_opt[long_opt_idx].has_arg = required_argument;
    long_opt[long_opt_idx].flag = NULL;
    long_opt[long_opt_idx++].val = 'd';

    long_opt[long_opt_idx].name = "accn";
    long_opt[long_opt_idx].has_arg = required_argument;
    long_opt[long_opt_idx].flag = NULL;
    long_opt[long_opt_idx++].val = 'a';

    long_opt[long_opt_idx].name = "jerk";
    long_opt[long_opt_idx].has_arg = required_argument;
    long_opt[long_opt_idx].flag = NULL;
    long_opt[long_opt_idx++].val = 'j';


    long_opt[long_opt_idx].name = "tune-accn";
    long_opt[long_opt_idx].has_arg = no_argument;
    long_opt[long_opt_idx].flag = NULL;
    long_opt[long_opt_idx++].val = 'A';

    long_opt[long_opt_idx].name = "tune-jerk";
    long_opt[long_opt_idx].has_arg = no_argument;
    long_opt[long_opt_idx].flag = NULL;
    long_opt[long_opt_idx++].val = 'J';

    long_opt[long_opt_idx].name = "nbins";
    long_opt[long_opt_idx].has_arg = required_argument;
    long_opt[long_opt_idx].flag = NULL;
    long_opt[long_opt_idx++].val = 'n';

    long_opt[long_opt_idx].name = "nsub";
    long_opt[long_opt_idx].has_arg = required_argument;
    long_opt[long_opt_idx].flag = NULL;
    long_opt[long_opt_idx++].val = 's';

    long_opt[long_opt_idx].name = "period-range";
    long_opt[long_opt_idx].has_arg = required_argument;
    long_opt[long_opt_idx].flag = &opt_flag;
    long_opt[long_opt_idx++].val = 1;

    long_opt[long_opt_idx].name = "period-step";
    long_opt[long_opt_idx].has_arg = required_argument;
    long_opt[long_opt_idx].flag = &opt_flag;
    long_opt[long_opt_idx++].val = 2;

    long_opt[long_opt_idx].name = "dm-range";
    long_opt[long_opt_idx].has_arg = required_argument;
    long_opt[long_opt_idx].flag = &opt_flag;
    long_opt[long_opt_idx++].val = 3;

    long_opt[long_opt_idx].name = "dm-step";
    long_opt[long_opt_idx].has_arg = required_argument;
    long_opt[long_opt_idx].flag = &opt_flag;
    long_opt[long_opt_idx++].val = 4;

    long_opt[long_opt_idx].name = "accn-range";
    long_opt[long_opt_idx].has_arg = required_argument;
    long_opt[long_opt_idx].flag = &opt_flag;
    long_opt[long_opt_idx++].val = 5;

    long_opt[long_opt_idx].name = "accn-step";
    long_opt[long_opt_idx].has_arg = required_argument;
    long_opt[long_opt_idx].flag = &opt_flag;
    long_opt[long_opt_idx++].val = 6;

    long_opt[long_opt_idx].name = "jerk-range";
    long_opt[long_opt_idx].has_arg = required_argument;
    long_opt[long_opt_idx].flag = &opt_flag;
    long_opt[long_opt_idx++].val = 7;

    long_opt[long_opt_idx].name = "jerk-step";
    long_opt[long_opt_idx].has_arg = required_argument;
    long_opt[long_opt_idx].flag = &opt_flag;
    long_opt[long_opt_idx++].val = 8;

    long_opt[long_opt_idx].name = "rms-size";
    long_opt[long_opt_idx].has_arg = required_argument;
    long_opt[long_opt_idx].flag = &opt_flag;
    long_opt[long_opt_idx++].val = 9;

    long_opt[long_opt_idx].name = "ascii-dump";
    long_opt[long_opt_idx].has_arg = no_argument;
    long_opt[long_opt_idx].flag = &opt_flag;
    long_opt[long_opt_idx++].val = 256;

    long_opt[long_opt_idx].name = "input-phcx";
    long_opt[long_opt_idx].has_arg = required_argument;
    long_opt[long_opt_idx].flag = &opt_flag;
    long_opt[long_opt_idx++].val = 257;

    long_opt[long_opt_idx].name = "output-phcx";
    long_opt[long_opt_idx].has_arg = required_argument;
    long_opt[long_opt_idx].flag = &opt_flag;
    long_opt[long_opt_idx++].val = 258;

    long_opt[long_opt_idx].name = "input-rawcand";
    long_opt[long_opt_idx].has_arg = required_argument;
    long_opt[long_opt_idx].flag = &opt_flag;
    long_opt[long_opt_idx++].val = 259;

    long_opt[long_opt_idx].name = "output-rawcand";
    long_opt[long_opt_idx].has_arg = required_argument;
    long_opt[long_opt_idx].flag = &opt_flag;
    long_opt[long_opt_idx++].val = 260;

    long_opt[long_opt_idx].name = 0;
    long_opt[long_opt_idx].has_arg = 0;
    long_opt[long_opt_idx].flag = 0;
    long_opt[long_opt_idx++].val = 0;

    pch_tune_state* state = (pch_tune_state*) malloc(sizeof (pch_tune_state));

    state->period = 1;
    state->acc = 0;
    state->jerk = 0;
    state->dm = 0;
    state->profile_nbins = 128;
    state->nsubints = 128;

    state->use_acc = 0;
    state->use_jerk = 0;

    state->period_range = NOTSETVAL;
    state->period_step = NOTSETVAL;
    state->dm_step = NOTSETVAL;
    state->dm_range = NOTSETVAL;
    state->acc_step = NOTSETVAL;
    state->acc_range = NOTSETVAL;
    state->jerk_step = NOTSETVAL;
    state->jerk_range = NOTSETVAL;

    state->rms_size = 10000;

    ascii_dump = 0;

    while ((c = getopt_long(argc, argv, args, long_opt, &long_opt_idx)) != -1) {
        // handle args
        switch (opt_flag) {
            case 1: // period-range
                state->period_range = atof(optarg) / 1000.0;
                opt_flag = 0;
                break;
            case 2: // period-step
                state->period_step = atof(optarg) / 1000.0;
                opt_flag = 0;
                break;
            case 3: // dm-range
                state->dm_range = atof(optarg);
                opt_flag = 0;
                break;
            case 4: // dm-step
                state->dm_step = atof(optarg);
                opt_flag = 0;
                break;
            case 5: // accn-range
                state->acc_range = atof(optarg);
                state->use_acc = 1;
                opt_flag = 0;
                break;
            case 6: // accn-step
                state->acc_step = atof(optarg);
                state->use_acc = 1;
                opt_flag = 0;
                break;
            case 7: // jerk-range
                state->jerk_range = atof(optarg);
                state->use_jerk = 1;
                opt_flag = 0;
                break;
            case 8: // jerk-step
                state->jerk_step = atof(optarg);
                state->use_jerk = 1;

                opt_flag = 0;
                break;
            case 9: // rms-size
                state->rms_size = atoi(optarg);
                opt_flag = 0;
                break;
            case 256:
                ascii_dump = 1;
                opt_flag = 0;
                break;
            case 257:
                if (have_input_rawcand == 1) {
                    fprintf(stderr, "Warning: Cannot input both phcx and rawcand, ignoring last\n");
                } else {
                    convert_from_bary = 0;
                    have_input_phcx = 1;
                    strcpy(input_phcx_filename, optarg);
                    input_phcx = read_phcx(input_phcx_filename);
                    pch_tune_copy_phcx_to_state(state, input_phcx);
                }
                opt_flag = 0;
                break;
            case 258:
                have_output_phcx = 1;
                strcpy(output_phcx_filename, optarg);
                opt_flag = 0;
                break;
            case 259:
                if (have_input_phcx == 1) {
                    fprintf(stderr, "Warning: Cannot input both phcx and rawcand, ignoring last\n");
                } else {
                    have_input_rawcand = 1;
                    convert_from_bary = 0;
                    strcpy(input_rawcand_filename, optarg);
                    input_phcx = read_phcx(input_rawcand_filename);
                    pch_tune_copy_phcx_to_state(state, input_phcx);
                }
                opt_flag = 0;
                break;
            case 260:
                have_output_rawcand = 1;
                strcpy(output_rawcand_filename, optarg);
                opt_flag = 0;
                break;
            default:
                opt_flag = 0;
                switch (c) {
                    case 'h':
                        pch_tune_print_usage();
                        exit(0);
                        break;
                    case 'p':
                        state->period = atof(optarg) / 1000.0;
                        break;
                    case 'P':
                        state->period = atof(optarg) / 1000.0;
                        // convert to bary...;
                        convert_from_bary = 1;
                        break;
                    case 'd':
                        state->dm = atof(optarg);
                        dm_set = 1;
                        break;
                    case 'a':
                        state->acc = atof(optarg);
                        break;
                    case 'j':
                        state->jerk = atof(optarg);
                        break;
                    case 'A':
                        state->use_acc = 1;
                        break;
                    case 'J':
                        state->use_jerk = 1;
                        break;
                    case 'n':
                        state->profile_nbins = atoi(optarg);
                        break;
                    case 's':
                        state->nsubints = atoi(optarg);
                        break;
                    case '?':
                        fprintf(stderr, "Error: Unknown option passed... exiting\n");
                        exit(1);
                        break;
                }
                break;
        }
    }

    if (optind >= argc) {
        printf("Error: Need a psrxml file to work on... (use --help for options)\n");
        exit(1);
    }

    readPsrXml(header, argv[optind]);

    state->header = header;
    state->nchans = header->numberOfChannels;
    double doppler_fac = 1;
    int telid = get_telid(header->telescope.pulsarhunterCode);
    if (telid == 0) {
        fprintf(stderr, "Warning: Unknown telescope used, cannot barycentre\n");
    } else {
        double epoch = ((double) header->mjdObs) / 86400.0 + (double) (header->timeToFirstSample) / 1e9;
        doppler_fac = barycentre_doppler_factor(epoch, telid, header->startCoordinate.ra, header->startCoordinate.dec);
    }
    if (convert_from_bary) {
        state->period *= doppler_fac; // bary->topo

    }


    if (state->period_range > NOTSETLIMIT) {
        state->period_range = state->period / 1000;
    }
    if (state->period_step > NOTSETLIMIT) {
        state->period_step = state->period_range / 100;
    }

    if (header->numberOfChannels == 1) {
        state->dm_range = 0;
        state->dm_step = 1;
    }else {
        if (state->dm_range > NOTSETLIMIT) {
            state->dm_range = 25;
        }
        if (state->dm_step > NOTSETLIMIT) {
            state->dm_step = 5;
        }
    }

    if (state->use_acc) {
        if (state->acc_range > NOTSETLIMIT) {
            state->acc_range = 100;
        }
        if (state->acc_step > NOTSETLIMIT) {
            state->acc_step = state->acc_range / 10;
        }
    }

    if (state->use_jerk) {
        if (state->jerk_range > NOTSETLIMIT) {
            state->jerk_range = 1;
        }
        if (state->jerk_step > NOTSETLIMIT) {
            state->jerk_step = state->jerk_range / 10;
        }
    }

    state->file_num = 0;
    state->data_file_name = header->files[0]->filename;
    if (state->profile_nbins > (int) (state->period
            / header->currentSampleInterval + 0.5)) {
        state->profile_nbins = (int) (state->period / header->currentSampleInterval
                + 0.5);
    }

    if (state->nsubints * state->period > header->actualObsTime) {
        fprintf(stderr, "Too many sub-ints for this period and tobs\n");
        exit(2);
    }

    printf("We have %d channels, %d subints and %d profile bins\n",
            state->nchans, state->nsubints, state->profile_nbins);

    /*
     * Do the folding of the sub-ints
     * 
     * 
     */

    pch_tune_scrunched_t* scrunched = pch_tune_make_scrunched(state);

    //	printf("RMS: %f\n",scrunched->rms);
    /*	for(int ch = 0; ch < state->nsubints; ch++){
     for(int bin = 0; bin < state->nbins; bin++){
     printf("%d %d %f\n",ch,bin,scrunched->dataCube[0][ch][bin]);
     }
     printf("\n");
     }*/

    /*
     * do the optimisation...
     * 
     */

    printf("Optimising period %12.9lf -> %12.9lf (%lf)\n", 1000 * (state->period
            - state->period_range), 1000.0 * (state->period + state->period_range),
            1000.0 * state->period_step);
    if (state->dm_range > 0)printf("Optimising dm %lf -> %lf (%lf)\n", (state->dm
            - state->dm_range), (state->dm + state->dm_range),
            state->dm_step);
    if (state->use_acc)printf("Optimising acc %lf -> %lf (%lf)\n", (state->acc
            - state->acc_range), (state->acc + state->acc_range),
            state->acc_step);
    if (state->use_jerk)printf("Optimising jerk %lf -> %lf (%lf)\n", (state->jerk
            - state->jerk_range), (state->jerk + state->jerk_range),
            state->jerk_step);

    pch_tune_optimise_result_t* opt_result = pch_tune_optimise_scrunched(state,
            scrunched);
    pch_tune_stack_result_t* sr =
            opt_result->stack_results[opt_result->best_dm_step];
    printf("Best Topo Period: %16.12lf ms\n", sr->best_period * 1000.0);
    printf("Best DM: %lf\n", opt_result->best_dm);

    printf("Best SNR: %lf\n", sr->best_snr);


    if (ascii_dump) {
        printf("Writing ascii dump file... ");
        fflush(stdout);

        // write out some helpful files for plotting
        FILE* ascii_file;
        ascii_file = fopen("pch_subints.asc", "w");
        for (int sub = 0; sub < state->nsubints; sub++) {
            for (int bin = 0; bin < state->profile_nbins; bin++) {
                fprintf(ascii_file, "%d %d %f\n", sub, bin,
                        opt_result->best_subints[sub][bin]);
            }
            fprintf(ascii_file, "\n");
        }
        fclose(ascii_file);

        ascii_file = fopen("pch_subbands.asc", "w");
        for (int ch = 0; ch < state->nchans; ch++) {
            float max = -1e16;
            float min = 1e16;
            for (int bin = 0; bin < state->profile_nbins; bin++) {
                float v = opt_result->best_subbands[ch][bin];
                if (v > max)max = v;
                if (v < min)min = v;
            }
            for (int bin = 0; bin < state->profile_nbins; bin++) {

                fprintf(ascii_file, "%d %d %f\n", ch, bin,
                        (opt_result->orig_subbands[ch][bin] - min) / (max - min));
            }
            fprintf(ascii_file, "\n");

        }
        fclose(ascii_file);

        ascii_file = fopen("pch_profile.asc", "w");
        for (int bin = 0; bin < state->profile_nbins; bin++) {
            fprintf(ascii_file, "%d %f\n", bin, sr->best_profile[bin]);
        }
        fclose(ascii_file);
        printf("Done\n");
    }

    if (have_output_phcx) {
        printf("Writing phcx file... ");
        fflush(stdout);
        if (input_phcx == NULL) {
            input_phcx = (phcx*) malloc(sizeof (phcx));

            input_phcx->header.centreFreq = header->centreFreqCh1 + header->freqOffset * header->numberOfChannels / 2;
            input_phcx->header.bandwidth = header->freqOffset * header->numberOfChannels;
            input_phcx->header.mjdStart = header->mjdObs;
            input_phcx->header.observationLength = header->actualObsTime;
            strcpy(input_phcx->header.sourceID, header->sourceName);
            strcpy(input_phcx->header.telescope, header->telescope.pulsarhunterCode);
            input_phcx->header.ra = header->startCoordinate.ra;
            input_phcx->header.dec = header->startCoordinate.dec;



            input_phcx->nsections = 1;
            input_phcx->sections = (phcx_section*) malloc(sizeof (phcx_section));
            input_phcx->sections->name = (char*) malloc(80);
            strcpy(input_phcx->sections->name, "UserDefined");
            input_phcx->sections->bestTopoPeriod = state->period;
            input_phcx->sections->bestBaryPeriod = state->period / doppler_fac;

            input_phcx->sections->bestDm = state->dm;
            input_phcx->sections->bestAccn = state->acc;
            input_phcx->sections->bestJerk = state->jerk;
            input_phcx->sections->tsamp = header->currentSampleInterval * 1e6;


        }

        if (input_phcx->sections->subbands == NULL) {
            input_phcx->sections->nsubbands = state->nchans;
            input_phcx->sections->nbins = state->profile_nbins;
            input_phcx->sections->subbands = opt_result->orig_subbands;
        }
        if (input_phcx->sections->subbands == NULL) {
            input_phcx->sections->nsubints = state->nsubints;
            input_phcx->sections->nbins = state->profile_nbins;
            input_phcx->sections->subbands = opt_result->orig_subints;
        }
        if (input_phcx->sections->subbands == NULL) {
            input_phcx->sections->nbins = state->profile_nbins;
            input_phcx->sections->pulseProfile = opt_result->orig_profile;
        }



        int i = input_phcx->nsections; // last section number.
        input_phcx->nsections++;
        input_phcx->sections = (phcx_section*) realloc(input_phcx->sections, sizeof (phcx_section) * input_phcx->nsections);
        input_phcx->sections[i].name = (char*) malloc(strlen(input_phcx->sections[i - 1].name) + 20);
        sprintf(input_phcx->sections[i].name, "%s-Tuned", input_phcx->sections[i - 1].name);
        input_phcx->sections[i].bestTopoPeriod = opt_result->stack_results[opt_result->best_dm_step]->best_period;
        input_phcx->sections[i].bestBaryPeriod = opt_result->stack_results[opt_result->best_dm_step]->best_period / doppler_fac;

        input_phcx->sections[i].bestAccn = opt_result->stack_results[opt_result->best_dm_step]->best_accn;
        input_phcx->sections[i].bestJerk = opt_result->stack_results[opt_result->best_dm_step]->best_jerk;
        input_phcx->sections[i].bestDm = opt_result->best_dm;
        input_phcx->sections[i].bestSnr = opt_result->stack_results[opt_result->best_dm_step]->best_snr;
        input_phcx->sections[i].bestWidth = opt_result->stack_results[opt_result->best_dm_step]->best_width;
        input_phcx->sections[i].tsamp = header->currentSampleInterval * 1e6;


        input_phcx->sections[i].nbins = state->profile_nbins;
        input_phcx->sections[i].nsubints = state->nsubints;
        input_phcx->sections[i].nsubbands = state->nchans;

        input_phcx->sections[i].pulseProfile = opt_result->stack_results[opt_result->best_dm_step]->best_profile;
        input_phcx->sections[i].subints = opt_result->best_subints;
        if (state->nchans > 1)input_phcx->sections[i].subbands = opt_result->best_subbands;

        memcpy(&(input_phcx->sections[i].snrBlock), opt_result->snr_block, sizeof (phcx_SNRBlock));


        write_phcx(output_phcx_filename, input_phcx);

        printf("Done\n");
    }

    // Free the memory!
    for (int ch = 0; ch < state->nchans; ch++) {
        for (int s = 0; s < state->nsubints; s++) {
            free(scrunched->dataCube[ch][s]);
        }
        free(scrunched->dataCube[ch]);
    }
    free(scrunched->dataCube);
    free(scrunched);

    for (int dm = 0; dm < opt_result->ndms; dm++) {
        free(opt_result->stack_results[dm]->best_profile);
        free(opt_result->stack_results[dm]);
    }

    for (int s = 0; s < state->nsubints; s++) {
        free(opt_result->best_subints[s]);
    }
    free(opt_result->best_subints);
    for (int ch = 0; ch < state->nchans; ch++) {
        free(opt_result->best_subbands[ch]);
    }
    free(opt_result->best_subbands);

    if (opt_result->snr_block != NULL) {
        phcx_snr_block* snr_block = opt_result->snr_block;
        free(snr_block->periodIndex);
        free(snr_block->dmIndex);
        free(snr_block->jerkIndex);
        free(snr_block->accnIndex);
        for (int a = 0; a < snr_block->ndm; a++) {
            for (int b = 0; b < snr_block->nperiod; b++) {
                for (int c = 0; c < snr_block->naccn; c++) {
                    free(snr_block->block[a][b][c]);
                }
                free(snr_block->block[a][b]);
            }
            free(snr_block->block[a]);
        }
        free(snr_block->block);
        free(snr_block);
    }
    free(opt_result);

    //	freePsrXml(header);
    free(state);


    printf("pch-tune complete\n");
}

void pch_tune_copy_phcx_to_state(pch_tune_state_t* state, phcx* infile) {

    state->period = infile->sections[infile->nsections - 1].bestTopoPeriod;
    state->dm = infile->sections[infile->nsections - 1].bestDm;
    state->acc = infile->sections[infile->nsections - 1].bestAccn;
    state->jerk = infile->sections[infile->nsections - 1].bestJerk;
    printf("Read candidate file: p %lf, d %lf, a %lf, j %lf\n", state->period
            * 1000.0, state->dm, state->acc, state->jerk);

}
