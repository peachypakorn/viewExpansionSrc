#ifndef __S4_FINALRECON__
#define __S4_FINALRECON__

#include "pyrrecon.h"

void S4_finalrecon(hls::stream< t_recon_complex>  &DFTOut,
		hls::stream< PIXEL_RAW>  &sigOut
		,hls::stream<ifft_out_t> &outf);



#endif
