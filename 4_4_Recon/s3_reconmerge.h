#ifndef __S3__RECON
#define __S3__RECON
#include "pyrrecon.h"

void S3_ReconMerge(hls::stream< t_nufft_output_complex > &H,
					hls::stream< t_nufft_output_complex > &L0,
				hls::stream< t_nufft_output_complex > &LA,
				hls::stream< t_nufft_output_complex > &LP,
		        hls::stream< t_recon_complex>  &imDFTOut);

#endif
