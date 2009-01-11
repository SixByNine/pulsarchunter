#include "pch-tune.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#define KDM 4.148808e-3

pch_tune_scrunched* pch_tune_make_scrunched(pch_tune_state_t* state) {
    // variable decliarations
    psrxml* header;
    dataFile* file;
    int i, j, bin, tmp = 0;
    char swap_chans;
    float*** profiles;

    header = state->header;
    file = header->files[state->file_num];
    swap_chans = header->freqOffset > 0;

    int samp_per_subint = header->numberOfSamples / state->nsubints;
    int samples_per_byte = 8 / file->bitsPerSample;
    int bytes_per_subint = samp_per_subint / samples_per_byte;
    //
    //	if (bytes_per_subint % file->blockLength) {
    //		bytes_per_subint = (int)(bytes_per_subint/file->blockLength)
    //				*file->blockLength;
    //		samp_per_subint = samples_per_byte * bytes_per_subint;
    //	}

    double time_per_subint = samp_per_subint * header->currentSampleInterval;

    const int profile_nbins = state->profile_nbins;
    const int nchan = state->nchans;
    const double sampleRate = header->currentSampleInterval;

    float** chans = (float**) malloc(sizeof (float*) * nchan);

    profiles = (float***) malloc(sizeof (float**) * nchan);
    for (i = 0; i < nchan; i++) {
        profiles[i] = (float**) malloc(sizeof (float*) * state->nsubints);
        for (j = 0; j < state->nsubints; j++) {
            profiles[i][j] = (float*) malloc(sizeof (float) * state->profile_nbins);
            for (bin = 0; bin < profile_nbins; bin++) {
                profiles[i][j][bin] = 0;
            }
        }
    }
    // prepare the data file...
    readPsrXMLPrepDataFile(file, state->data_file_name);

    double epoch = header->actualObsTime / 2.0;
    epoch += header->actualObsTime / 2.0; // MJD at centre of observation.

    uint64 start = 0;
    uint64 end = samp_per_subint;

    int rms_bins = 0;
    double rms = 0;
    double sum = 0;
    double* rms_ch = (double*) malloc(sizeof (double) * nchan);
    double* sum_ch = (double*) malloc(sizeof (double) * nchan);
    for (int j = 0; j < nchan; j++) {
        rms_ch[j] = 0;
        sum_ch[j] = 0;
    }
    const int rms_limit = state->rms_size;

    int excess_bytes_length = 0;
    char* excess_bytes = (char*) malloc(sizeof (char) * file->blockLength);
    int read = 0;

    for (i = 0; i < state->nsubints; i++) {
        fprintf(stdout, "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
        fprintf(stdout, "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
        fprintf(stdout, "Folding sub-integration %d/%d", i + 1, state->nsubints);
        fflush(stdout);

        /*
         * We read some data into an array.
         * 
         * 
         */
        float* sub_data = (float*) malloc(sizeof (float) * samp_per_subint
                * header->numberOfChannels);
        unsigned char* sub_bytes = (unsigned char*) malloc(sizeof (unsigned char)
                * samp_per_subint * header->numberOfChannels);

        int nbytes_to_read = samp_per_subint * file->bitsPerSample
                * header->numberOfChannels / 8;
        int blocks_read = 0;
        int subint_bytes_read = 0;
        if (excess_bytes_length > 0) {
            int sz = excess_bytes_length;
            if (sz > nbytes_to_read)
                sz = nbytes_to_read;
            memcpy(sub_bytes, excess_bytes, sz);
            nbytes_to_read -= sz;
            excess_bytes_length = 0;

        }
        while (nbytes_to_read > file->blockLength) {
            read = readPsrXmlNextDataBlockIntoExistingArray(file, sub_bytes
                    + file->blockLength * blocks_read);
            if (read < 0) {
                // we have reached the end of the file...
                fprintf(stderr, "\n\nReached the end of file...\n");

                nbytes_to_read = 0;
                samp_per_subint = (int) (subint_bytes_read * samples_per_byte
                        / header->numberOfChannels);
                break;
            }
            nbytes_to_read -= read;
            subint_bytes_read += read;
            blocks_read++;
        }
        //		fprintf(stderr,"Read %d bytes\n", subint_bytes_read);
        // we now may have a few bytes to fill our array, but less than one block worth.
        if (nbytes_to_read > 0) {
            unsigned char* block;
            readPsrXmlNextDataBlock(file, &block);

            memcpy(sub_bytes + file->blockLength*blocks_read, block,
                    nbytes_to_read);
            excess_bytes_length = file->blockLength - nbytes_to_read;
            memcpy(excess_bytes, block + nbytes_to_read, excess_bytes_length);
            nbytes_to_read = 0;
        }
        // this converts to the raw data to nicely arranged floats
        unpackDataChunk(sub_bytes, sub_data, header, state->file_num,
                samp_per_subint, start, end, swap_chans);
        unpackToChannels(sub_data, chans, nchan, samp_per_subint);

        // free any intermediate arrays...
        free(sub_bytes);

        if (rms_bins < rms_limit) {
            // we want to add to the rms sum
            for (bin = 0; rms_bins < rms_limit && bin < samp_per_subint; bin++) {
                for (j = 0; j < nchan; j++) {
                    float v = chans[j][bin];
                    rms_ch[j] += v*v;
                    sum_ch[j] += v;
                }
                rms_bins++;
            }
        }

        /*
         * Now we have a nice float array, we should start folding...
         * 
         * 
         */
        double toff_ref = i * time_per_subint + epoch;

        int** count = (int**) malloc(sizeof (int*) * nchan);
        for (j = 0; j < nchan; j++) {
            count[j] = (int*) malloc(sizeof (int) * state->profile_nbins);
            for (bin = 0; bin < profile_nbins; bin++) {
                count[j][bin] = 0;

            }
        }
        /*
         * Commented section writes out in 'reader' style to check file is read ok.
         * 
         */
        //		printf("\n");
        //		for (bin = 0; bin < samp_per_subint && bin < 1000; bin++) {
        //			printf("%f", tmp*header->currentSampleInterval);
        //			tmp++;
        //			for (j = 0; j < nchan; j++) {
        //				printf(" %f", chans[j][bin]);
        //			}
        //			printf("\n");
        //		}

        double target_bin_d = 0;
        int target_bin;
        for (j = 0; j < nchan; j++) {
            double natural_time_passed = toff_ref;
            double pulsar_time_passed;
            for (bin = 0; bin < samp_per_subint; bin++) {

                // This addjusts the 'flow of time' to correct for accn/jerk
                pulsar_time_passed = pch_tune_get_adjusted_time(state,
                        natural_time_passed);

                // target_bin is the bin into which to write the data. _d is a floating point representation.
                target_bin_d = (((pulsar_time_passed) / state->period))
                        * profile_nbins + 0.5;

                target_bin = (((int) target_bin_d) % profile_nbins);

                // correct if we have gone below 0.
                if (target_bin_d < 0) {
                    target_bin--;
                    target_bin += profile_nbins;
                }

                profiles[j][i][target_bin] += chans[j][bin];
                count[j][target_bin] += 1;
                natural_time_passed += sampleRate;
            }
        }

        double countoff = (double) samp_per_subint / (double) profile_nbins;
        for (j = 0; j < nchan; j++)
            for (bin = 0; bin < profile_nbins; bin++) {
                //printf("%f %f %d\n",profiles[j][i][bin],countoff,count[j][bin]);
                profiles[j][i][bin] *= (countoff / (float) count[j][bin]);
            }

        // free arrays associated with this subint...
        free(sub_data);
    }

    printf("\nCreated all sub-integrations.\n");

    // free remaining arrays
    free(excess_bytes);
    free(chans);

    for (int j = 0; j < nchan; j++) {
        rms_ch[j] /= rms_bins;
        sum_ch[j] /= rms_bins;

        double meansq = sum_ch[j] * sum_ch[j];

        rms_ch[j] = rms_ch[j] - meansq;
                rms += rms_ch[j];
    }
    rms = sqrt(rms);

    printf("RMS: %lf (from %d bins)\n", rms, rms_bins);
    pch_tune_scrunched_t* result =
            (pch_tune_scrunched_t*) malloc(sizeof (pch_tune_scrunched_t));
    result->dataCube = profiles;
    result->rms = rms;
    return result;

}

pch_tune_optimise_result_t* pch_tune_optimise_scrunched(
        pch_tune_state_t* state, pch_tune_scrunched_t* scrunched) {

    psrxml* header = state->header;

    const int nchan = state->nchans;

    double binWidth = state->period / (double) state->profile_nbins;

    //	double numFoldsPerSubint = (header->numberOfSamples/state->period)
    //			/state->nsubints;
    //	double npulsebin = state->period/header->currentSampleInterval;

    double numFoldsPerSubint = (header->actualObsTime / state->period)
            / state->nsubints;
    double npulsebin = state->period / header->currentSampleInterval;

    double rmss = scrunched->rms * sqrt(numFoldsPerSubint * (npulsebin
            / (double) state->profile_nbins)); // rms in a subint


    double flo = header->centreFreqCh1 / 1000.0; //GHz
    double fhi = (header->centreFreqCh1 + header->numberOfChannels
            * header->freqOffset) / 1000.0; //GHz
    if (flo > fhi) {
        double fswap = fhi;
        fhi = flo;
        flo = fswap;
    }

    double minDdm = binWidth / (KDM * (1.0 / (flo * flo) - 1.0 / (fhi * fhi)));

    double deltaDm = state->dm_step;
    if (deltaDm < minDdm) {
        state->dm_range = ((minDdm / deltaDm) * state->dm_range);
        deltaDm = minDdm;
        if (nchan > 1)
            fprintf(stderr, "Dm Step below min threshold: reset to %lf\n"
                , deltaDm);
        state->dm_step = deltaDm;
    }

    double dmcoursness = (deltaDm * (KDM * (1.0 / (flo * flo) - 1.0 / (fhi * fhi)))) / (nchan * binWidth);

    int ndms = (int) (state->dm_range / deltaDm);
    if (nchan == 1) {
        ndms = 0;
        deltaDm = 0;
    }

    int dmStep = 0;
    int ptr = 0;
    fprintf(stderr, "Dedispersing subbands\n");
    for (int b = 0; b < nchan; b++) {

        double fch = (header->centreFreqCh1 + b * header->freqOffset) / 1000.0; //GHz		

        double timeOffset = state->dm * (KDM * (1.0 / (fch * fch) - 1.0
                / (fhi * fhi)));

        double coreDmOffset = (timeOffset / state->period) * state->profile_nbins;

        ptr = (int) ((coreDmOffset));

        /*
         *  this rotates the arrays to dedisperse them...
         */
        for (int i = 0; i < state->nsubints; i++) {

            float* tmpArr = (float*) malloc(sizeof (float) * state->profile_nbins);

            memcpy(tmpArr, scrunched->dataCube[b][i], sizeof (float)
                    * state->profile_nbins);
            for (int j = 0; j < state->profile_nbins; j++) {
                while (ptr >= state->profile_nbins)
                    ptr -= state->profile_nbins;
                while (ptr < 0)
                    ptr += state->profile_nbins;

                scrunched->dataCube[b][i][j] = tmpArr[ptr];

                ptr++;
            }
        }
    }

    /*
     * Now do some searching...
     * 
     * 
     */
    // this stores the results of the search...
    pch_tune_stack_result_t** results =
            (pch_tune_stack_result_t**) malloc(sizeof (pch_tune_stack_result_t*)
            *(2 * ndms + 1));
    pch_tune_stack_result_t* sr_best = NULL;
    float** bestSubints = NULL;
    float bestDm = 0;
    int bestDmStep = 0;
    phcx_snr_block* snr_block = NULL;
    for (int dms = -ndms; dms <= ndms; dms++) {

        double dmf = dms*dmcoursness;
        float** subints = pch_tune_dedisperse(state, dmf, scrunched->dataCube);

        pch_tune_stack_result_t* stack_slide_result = pch_tune_stack_slide(
                state, subints, &snr_block, rmss, dmStep);

        results[dmStep] = stack_slide_result;

        if (sr_best == NULL || sr_best->best_snr < stack_slide_result->best_snr) {
            sr_best = stack_slide_result;
            if (bestSubints != NULL)
                free(bestSubints); // if we had one already, free the memory
            bestSubints = subints;
            bestDm = dmf;
            bestDmStep = dmStep;
        } else {
            free(subints); // if this was the best, we want to keep the memory.
        }

        dmStep++;
    }
    printf("\nDone stack-slide\n");

    printf("Making optimised sub-bands\n");
    float** sbands = (float**) malloc(sizeof (float*) * nchan);
    float** old_sbands = (float**) malloc(sizeof (float*) * nchan);

    for (int ch = 0; ch < nchan; ch++) {
        old_sbands[ch] = (float*) malloc(sizeof (float) * state->profile_nbins);
        sbands[ch] = (float*) malloc(sizeof (float) * state->profile_nbins);
        for (int bin = 0; bin < state->profile_nbins; bin++) {
            sbands[ch][bin] = 0;
            old_sbands[ch][bin] = 0;
        }
    }
    for (int ch = 0; ch < nchan; ch++) {
        double b2 = (double) ch - nchan / 2.0;

        int ddptr = (int) ((bestDm) * b2);
        double s2 = -(state->nsubints / 2.0);
        const int nbin = state->profile_nbins;
        for (int s = 0; s < state->nsubints; s++) {
            s2 += 1;
            int pptr = (int) (sr_best->best_pfactor * s2 + sr_best->best_pdfactor
                    * s2 * s2 + sr_best->best_pddfactor * s2 * s2 * s2);
            int ptr = ddptr + pptr;
            if (ptr >= nbin)
                ptr -= nbin;
            while (ptr < 0)
                ptr += nbin;
            while (pptr < 0)
                pptr += nbin;

            for (int i = 0; i < nbin; i++) {

                while (ptr >= nbin)
                    ptr -= nbin;
                while (pptr >= nbin)
                    pptr -= nbin;
                sbands[ch][i] += scrunched->dataCube[ch][s][ptr];

                old_sbands[ch][i] += scrunched->dataCube[ch][s][pptr];

                ptr++;
                pptr++;
            }

        }
    }

    pch_tune_optimise_result_t
            * ret =
            (pch_tune_optimise_result_t*) malloc(sizeof (pch_tune_optimise_result_t));
    ret->stack_results = results;
    ret->ndms = (2 * ndms + 1);
    ret->best_dm_step = bestDmStep;
    // Now make bestDm have the best dm...
    if (nchan > 1) {
        double deltaT = bestDm * nchan * binWidth;
        deltaDm = deltaT / (KDM * (1.0 / (flo * flo) - 1.0 / (fhi * fhi)));
        ret->best_dm = state->dm + deltaDm;
    } else {
        ret->best_dm = state->dm;
    }
    ret->best_subints = bestSubints;
    ret->best_subbands = sbands;
    ret->orig_subbands = old_sbands;

    ret->snr_block = snr_block;
        printf("Optimisation Complete\n");

    return ret;
}

pch_tune_stack_result_t* pch_tune_stack_slide(pch_tune_state_t* state,
        float** subints, phcx_snr_block** snr_block_ptr, float rmss, int dmStep) {
    phcx_snr_block* snr_block = *snr_block_ptr;
    const int nsub = state->nsubints;
    const int nbins = state->profile_nbins;
    const float tobs = state->header->actualObsTime;
    pch_tune_stack_result_t* result =
            (pch_tune_stack_result_t*) malloc(sizeof (pch_tune_stack_result_t));
    result->best_profile = (float*) malloc(sizeof (float) * nbins);
    result->best_snr = -1000000;

    int nperiodTrials = (int) (state->period_range / state->period_step);

    int nAccnTrials = 0;

    if (state->use_acc)
        nAccnTrials = (int) (state->acc_range / state->acc_step);
    int nJerkTrials = 0;

    if (state->use_jerk)
        nJerkTrials = (int) (state->jerk_range / state->jerk_step);

    int nDmTrials = (int) (state->dm_range / state->dm_step);

    double binWidth = state->period / (double) nbins;

    int pTrial = 0;
    int accnTrial = 0;
    int jerkTrial = 0;

    double bestPFactor = 0;
    double bestPdotFactor = 0;
    double bestPddotFactor = 0;

    // rms in profile
    double rmsp = rmss * sqrt(nsub);
    double sHalf = (nsub / 2.0);


    double pcoarseness = (state->period_step * tobs) / (state->period * nsub
            * binWidth);

    double pdConversion = (tobs * tobs) / (2.0 * nsub * nsub * binWidth * SPEED_OF_LIGHT);
    double pddConversion = (tobs * tobs * tobs) / (6.0 * nsub * nsub * nsub * binWidth
            * SPEED_OF_LIGHT);


    int dm_step_trials = (1 + 2 * nperiodTrials)*(1 + 2 * nAccnTrials) *(1 + 2
            * nJerkTrials);
    int p_step_trials = (1 + 2 * nAccnTrials) *(1 + 2
            * nJerkTrials);
    int total_trials = dm_step_trials * (1 + 2 * nDmTrials);

    if (snr_block == NULL) {
        printf("Searched 0/%d trials", total_trials);
        fflush(stdout);
        *snr_block_ptr = (phcx_snr_block*) malloc(sizeof (phcx_snr_block));
        snr_block = *snr_block_ptr;

        snr_block->nperiod = 2 * nperiodTrials + 1;
        snr_block->ndm = 2 * nDmTrials + 1;
        snr_block->naccn = 2 * nAccnTrials + 1;
        snr_block->njerk = 2 * nJerkTrials + 1;

        snr_block->periodIndex = (double*) malloc(sizeof (double)
                *(snr_block->nperiod));
        int i = 0;
        for (int ps = -nperiodTrials; ps <= nperiodTrials; ps++) {
            double periodFactor = (double) ps * (double) pcoarseness;
            snr_block->periodIndex[i] = state->period + ps * state->period_step;
            i++;
        }
        snr_block->accnIndex = (double*) malloc(sizeof (double)
                *(snr_block->naccn));
        i = 0;
        for (int acs = -nAccnTrials; acs <= nAccnTrials; acs++) {
            double accnOffset = acs * state->acc_step;
            snr_block->accnIndex[i] = state->acc + accnOffset;
            i++;
        }

        snr_block->jerkIndex = (double*) malloc(sizeof (double)
                *(snr_block->njerk));
        i = 0;
        for (int jes = -nJerkTrials; jes <= nJerkTrials; jes++) {
            double jerkOffset = jes * state->jerk_step;
            snr_block->jerkIndex[i] = state->jerk + jerkOffset;
            i++;
        }

        snr_block->dmIndex = (double*) malloc(sizeof (double) *(snr_block->ndm));
        if (snr_block->ndm > 1) {
            i = 0;
            for (int dms = -nDmTrials; dms <= nDmTrials; dms++) {
                double dmOffset = dms * state->dm_step;
                snr_block->dmIndex[i] = state->dm + dmOffset;
                i++;
            }
        } else {
            snr_block->dmIndex[0] = state->dm;
        }


        snr_block->block = (double****) malloc(sizeof (double***)
                *(snr_block->ndm));
        for (int a = 0; a < snr_block->ndm; a++) {
            snr_block->block[a] = (double***) malloc(sizeof (double**)
                    *(snr_block->nperiod));
            for (int b = 0; b < snr_block->nperiod; b++) {
                snr_block->block[a][b] = (double**) malloc(sizeof (double*)
                        *(snr_block->naccn));
                for (int c = 0; c < snr_block->naccn; c++) {
                    snr_block->block[a][b][c] = (double*) malloc(sizeof (double)
                            * snr_block->njerk);
                    for (int d = 0; d < snr_block->njerk; d++) {
                        snr_block->block[a][b][c][d] = 0;
                    }
                }
            }
        }
    }
    // we have now made the snrblock arrays.


    // Search periods
    pTrial = -1;
    for (int ps = -nperiodTrials; ps <= nperiodTrials; ps++) {
        pTrial++;
        double periodFactor = (double) ps * (double) pcoarseness;
        
        // Search pdots
        accnTrial = -1;
        for (int acs = -nAccnTrials; acs <= nAccnTrials; acs++) {
            accnTrial++;

            double accnOffset = acs * state->acc_step;

            double pdotFactor = pdConversion*accnOffset;

            // Search pddots
            jerkTrial = -1;
            for (int jes = -nJerkTrials; jes <= nJerkTrials; jes++) {
                jerkTrial++;

                double jerkOffset = jes * state->jerk_step;

                double pddotFactor = pddConversion*jerkOffset;

                float* profile = (float*) malloc(sizeof (float) * nbins);
                for (int i = 0; i < nbins; i++)
                    profile[i] = 0;
                double s2 = -sHalf;

                for (int s = 0; s < nsub; s++) {
                    s2 += 1;
                    int ptr = (int) (periodFactor * s2 + pdotFactor * s2 * s2
                            + pddotFactor * s2 * s2 * s2);
                    if (ptr >= nbins)
                        ptr -= nbins;
                    while (ptr < 0)
                        ptr += nbins;

                    for (int i = 0; i < nbins; i++) {
                        while (ptr >= nbins)
                            ptr -= nbins;
                        profile[i] += subints[s][ptr];

                        ptr++;
                    }

                }
                // Now calculate the SNR for this trial...
                double snr = 0;
                double width = 0;
                pch_tune_smooth(profile, 32, rmsp, nbins, &snr, &width);
                snr_block->block[dmStep][pTrial][accnTrial][jerkTrial] = snr;
                if (snr > result->best_snr) {

                    result->best_snr = snr;
                    result->best_width = width;
                    result->best_period = state->period + ps * state->period_step;

                    result->best_accn = state->acc + accnOffset;
                    result->best_jerk = state->jerk + jerkOffset;

                    result->best_p_trial = pTrial;
                    result->best_acc_trial = accnTrial;
                    result->best_jerk_trial = jerkTrial;

                    result->best_pfactor = periodFactor;
                    result->best_pdfactor = pdotFactor;
                    result->best_pddfactor = pddotFactor;

                    memcpy(result->best_profile, profile, nbins * sizeof (float));

                }
                free(profile);
            }
        }
        if(state->use_acc){
              printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
    printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
    printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");

    printf("Searched %d/%d trials", (dmStep) * dm_step_trials + (pTrial+1)*p_step_trials, total_trials);
        }
    }

    printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
    printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
    printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");

    printf("Searched %d/%d trials", (dmStep + 1) * dm_step_trials, total_trials);

    fflush(stdout);
    return result;
}

void pch_tune_smooth(float* pr, int maxwidth, double rmsp, int nbin,
        double* snr, double* width) {

    int kwmax = 0;
    double snrmax = 0, smmax = 0;
    int ksm, kw, ja;
    double s, sn, smax;

    //  remove baseline
    ksm = (int) (nbin / 2.5 + 0.5);
    smax = 0;
    for (int k = 0; k < ksm; k++) {
        smax = smax + pr[k % nbin];
    }

    for (int j = 0; j < nbin; j++) {
        s = 0.0;
        for (int k = 0; k < ksm; k++) {
            s = s + pr[(j + k) % nbin];
        }
        if (s < smax)
            smax = s;
    }
    smax = smax / ksm;
    for (int j = 0; j < nbin; j++) {
        pr[j] = pr[j] - smax;
    }
    //--------------------------------------


    kw = 1;
    //for(int nn = 0; nn < 1000; nn++){
    while (kw <= maxwidth && kw < nbin / 2) {

        s = 0.0;

        for (int k = 0; k < kw; k++) {
            s += pr[k];
        }

        smax = s;
        // New convolve
        for (int j = 0; j < nbin; j++) {

            ja = (j + kw) % nbin;
            s -= pr[j];
            s += pr[ja];
            if (s > smax)
                smax = s;
        }
        sn = smax / (rmsp * sqrt(kw * (1.0 + (float) kw / (float) nbin)));

        //sn=smax/(rmsp*(kw + (float)(kw)/(float)nbin));
        if (sn > snrmax) {
            snrmax = sn;
            kwmax = kw;
            smmax = smax / kw;

        }
        // Double width
        kw *= 2;
    }
    *snr = snrmax;
    *width = (double) kwmax / (double) nbin;
}

float** pch_tune_dedisperse(pch_tune_state_t* state, double dm,
        float*** profiles) {
    const int nsub = state->nsubints;
    const int nbin = state->profile_nbins;
    const int nchan = state->nchans;
    float** subints = (float**) malloc(sizeof (float*) * nsub);
    for (int i = 0; i < nsub; i++) {
        subints[i] = (float*) malloc(sizeof (float) * nbin);
        for (int j = 0; j < nbin; j++) {
            subints[i][j] = 0;
        }
    }

    for (int b = 0; b < nchan; b++) {
        double b2 = (double) b - nchan / 2.0;

        int ptr = (int) ((dm) * b2);

        for (int i = 0; i < nsub; i++) {

            for (int j = 0; j < nbin; j++) {
                while (ptr >= nbin)
                    ptr -= nbin;
                while (ptr < 0)
                    ptr += nbin;

                subints[i][j] += profiles[b][i][ptr];
                ptr++;
            }
        }

    }
    return subints;
}

double pch_tune_get_adjusted_time(pch_tune_state_t* state, double natural_time) {
    return natural_time - state->acc * natural_time * natural_time / (2.0
            * SPEED_OF_LIGHT) - state->jerk * natural_time * natural_time
            * natural_time / (6.0 * SPEED_OF_LIGHT);
}
