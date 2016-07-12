#ifndef __FORWARDFFTALL__
#define __FORWARDFFTALL__

#include "ap_fixed.h"
#include "hls_fft.h"

#include "common_type.h"

#include <complex>

void s2_fft_all_stream(
		hls::stream< t_nufft_output_complex > &inH,
		hls::stream< t_nufft_output_complex > &inL0,
		hls::stream< t_nufft_output_complex > &inLA,
		hls::stream< t_nufft_output_complex > &outH,
		hls::stream< t_nufft_output_complex > &outL0,
		hls::stream< t_nufft_output_complex > &outLA
		);

#endif
