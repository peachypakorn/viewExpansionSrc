
#include "s4_finalrecon.h"


void S4_finalrecon(hls::stream< t_recon_complex>  &DFTOut,
		hls::stream< PIXEL_RAW>  &sigOut
		,hls::stream<ifft_out_t> &outf) {
#pragma HLS INTERFACE axis port=outf
#pragma HLS DATA_PACK variable=outf
#pragma HLS INTERFACE axis depth=512 port=DFTOut

	// Enable this if instantiated this alone
#pragma HLS INTERFACE axis depth=512 port=sigOut

#pragma HLS data_pack variable=DFTOut
#pragma HLS data_pack variable=sigOut

#pragma HLS DATAFLOW

	pyr_recon_ifft(DFTOut, sigOut,outf);

}
