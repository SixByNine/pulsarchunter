#include "pch-seek.h"
#include <string.h>

void pch_seek_init_operations(pch_seek_operations_t* operations){
	operations->dump_tim=0;
	operations->fft_input=0;
	operations->form_amplitudes=0;
	operations->twiddle_amplitudes=0;
	operations->dump_amplitudes=0;
	operations->dump_phases=0;
	operations->dump_normalised=0;
	operations->dump_harmfolds=0;

	operations->hist_tim=0;
	operations->hist_amplitudes=0;
	operations->hist_normalised=0;
	operations->hist_harmfolds=0;

	operations->phase_fit=0;
	operations->search_chans=0;
	operations->fscrunch=0;
	operations->normalise_median=0;
	operations->normalise_powerlaw=0;
	operations->normalise_agl=0;




	operations->harmfold_simple=0;
	operations->harmfold_smart=0;
	operations->search_amplitudes=0;
	operations->write_prd=0;
	operations->write_presto_fft=0;
	operations->recon_add=0;
	operations->recon_ralph=0;
	operations->giant_search = 0;
	operations->append_output=0;
	operations->use_sigproc_zapfile=0;

	operations->ndm=0;
	operations->dmtrials=NULL;
	operations->harmfolds=NULL;
	operations->nharms=0;
	strcpy(operations->prdfile,"out.prd");
	operations->amp_thresh=5;
	operations->hfold_bonus_factor=0.2;
	strcpy(operations->giantfile,"giant.sp");
	strcpy(operations->zapfile,"zap.file");

	strcpy(operations->presto_fft_file,"presto");
}

