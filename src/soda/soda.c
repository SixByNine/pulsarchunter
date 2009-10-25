#include <stdlib.h>
#include "soda_dedisp.h"
#include <psrxml.h>
#include <stdio.h>
#include <math.h>

void dedisperse_chunk(SODA_TYPE** in_data, SODA_TYPE** dedisp_data,int nblocks);
SODA_TYPE** repack_chunk(unsigned char* rawchunk, int nblock, psrxml* header, char swapChannels);
void write_chunk(SODA_TYPE** dedisp_data, int nblocks);
void output_repack(SODA_TYPE** dedisp_data, int nblocks);

int main(int argc, char** argv){
	psrxml* header;
	unsigned char* dat;
	SODA_TYPE **in_data, **dedisp_data;
	int i,c,d,nblock,nsamps,read,nblocks_to_read,pack_nblock,nbytes,totread, pad_blocks;
	dataFile* file;

	// some arbatry chunk size
	nblock=13000*1024;

	// Read the PSRXML header
	header = (psrxml*) malloc(sizeof (psrxml));
	readPsrXml(header, argv[1]);

	// Do some checks that the data are processable with this SODA instance
	if (header->numberOfChannels != SODA_NCHAN){
		fprintf(stderr,"ERROR, nchans must be %d for this SODA instance\n",SODA_NCHAN);
		exit(1);
	}

	if (sizeof(SODA_TYPE) != SODA_COMP){
		fprintf(stderr,"ERROR, SODA_TYPE has the wrong size\n");
		exit(1);
	}

	if (sizeof(SODA_PAD_TYPE) != SODA_PAD){
		fprintf(stderr,"ERROR, SODA_PAD_TYPE has the wrong size\n");
		exit(1);
	}


	// How many blocks are required to pad the end of each chunk
	pad_blocks=ceil(SODA_MAXSAMP/(float)SODA_BS);

	// prepair to read the file
	file=header->files[0];
	readPsrXMLPrepDataFile(file, file->filename);

	// This is the number of bytes we need to read from the file for one chunk
	nbytes=((nblock+pad_blocks)*SODA_BS)*header->numberOfChannels*file->bitsPerSample/8;

	// this is the number of input blocks to read
	nblocks_to_read=ceil(nbytes/(float)header->files[0]->blockLength);

	printf("READ %d\n",nblocks_to_read);
	// Make space in RAM.
	dat=(unsigned char*)malloc(sizeof(unsigned char*)*header->files[0]->blockLength*nblocks_to_read);

	// Read a chunk worth of data.
	totread=0;
	i=0;
	while(totread < nbytes){
		read = readPsrXmlNextDataBlockIntoExistingArray(file, dat+header->files[0]->blockLength*i);
		if(read < 1){
			break;
		}
		totread+=read;
		i++;
	}

	int nblocks_read = (int)(totread/(SODA_BS*header->numberOfChannels*file->bitsPerSample/8));
	if (nblocks_read < nblock-pad_blocks)nblock=nblocks_read-pad_blocks;
	printf("%f %d\n",totread/(float)(SODA_BS*header->numberOfChannels*file->bitsPerSample/8),nblock);
	
	// Re-pack the newly read data
	printf("REPACK (read %d bytes)\n",totread);
	nsamps = (totread * (8/file->bitsPerSample)) / header->numberOfChannels;
	in_data=repack_chunk(dat,nsamps,header,0);

	pack_nblock=nblock*SODA_PAD/SODA_COMP;

	// if the data read wasn't sufficient, we are probably at EOF and therefore shrink the number of blocks
	// used to process
	printf("DEDISP\n");
	// do the dedispersion!
	dedisp_data=(SODA_TYPE**)malloc(sizeof(SODA_TYPE*)*SODA_NDM);
	for (d = 0; d < SODA_NDM  ; d++){
		dedisp_data[d]=(SODA_TYPE*)malloc(sizeof(SODA_TYPE)*SODA_BS*pack_nblock);
	}
	// here we pretend we have less blocks because we have packed many samples into one SODA_TYPE.
	dedisperse_chunk(in_data,dedisp_data,pack_nblock);

	// Write out the data.
	// @TODO!


	printf("PACK FOR OUTPUT\n");
	output_repack(dedisp_data,nblock);
	printf("WRITE\n");
	write_chunk(dedisp_data,nblock);

	printf("DONE\n");
	free(dat);
}


void output_repack(SODA_TYPE** dedisp_data, int nblocks){
	int nbin;
	int c,b;
	int factor,off;
	char* output;
	SODA_PAD_TYPE* input;
	factor=16;
	off=-128;
	nbin=SODA_BS*nblocks;
	for(c=0; c < SODA_NDM; c++){
		output = (char*) (dedisp_data[c]);
		input = (SODA_PAD_TYPE*) (dedisp_data[c]);
		for (b=0; b< nbin; b++){
			output[b]=(char)(input[b]/factor + off);
		}
	}
}


void write_chunk(SODA_TYPE** dedisp_data, int nblocks){
	FILE *file;
	char str[256];
	int d;
	for (d=0; d < SODA_NDM; d++){
		sprintf(str,"%d.raw",d);
		file = fopen(str,"w");
		fwrite(dedisp_data[d],1, nblocks*SODA_BS, file);
		fclose(file);
	}
}

void dedisperse_chunk(SODA_TYPE** in_data, SODA_TYPE** dedisp_data,int nblocks){
	int i,c,d;
	SODA_TYPE** inarr;
	SODA_TYPE** outarr;
	for (i = 0; i < nblocks; i++){
		inarr=(SODA_TYPE**) malloc(sizeof(SODA_TYPE*)*SODA_NCHAN);
		for (c = 0; c < SODA_NCHAN; c++){
			inarr[c]=in_data[c]+i*SODA_BS;
		}	
		outarr=(SODA_TYPE**) malloc(sizeof(SODA_TYPE*)*SODA_NDM);
		for (d = 0; d < SODA_NDM ; d++){
			outarr[d]=dedisp_data[d]+i*SODA_BS;
		}	
		
		soda_dedisp(outarr,inarr);
		free(inarr);
		free(outarr);
	}

}

SODA_TYPE** repack_chunk(unsigned char* rawchunk, int nsamp, psrxml* header, char swapChannels){
	int c;
	int fileNum=0;
	unsigned short* outshort=(unsigned short*)malloc(sizeof(short) * nsamp*SODA_NCHAN);
	unpackDataChunk_1to8bit_toshort(rawchunk, outshort,
			header->files[fileNum]->bitsPerSample,
			header->numberOfChannels, nsamp,
			header->files[fileNum]->firstSampleIsMostSignificantBit,
			header->files[fileNum]->isSigned,
			header->files[fileNum]->isChannelInterleaved, 0,
			nsamp, swapChannels);

	SODA_TYPE** out = (SODA_TYPE**)malloc(sizeof(SODA_TYPE*)*SODA_NCHAN);
	for (c=0; c < SODA_NCHAN; c++){
		out[c] = (SODA_TYPE*) &(outshort[nsamp*c]);
	}

	return out;
}


