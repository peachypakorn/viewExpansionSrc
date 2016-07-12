#include "nufft.h"
#include "nufft_cosfilter.h"




void nufft_cosfilter1(hls::stream< t_nufft_output_complex > nufftIn[1],
					 hls::stream< t_nufft_output_complex > nufftOut[1],
					 const int nL, const int m) {

#pragma HLS data_pack variable=nufftIn
#pragma HLS data_pack variable=nufftOut

#pragma HLS DATAFLOW
	//nufft_cosfilter<1>(nufftIn, nufftOut, nL, m);
	nufft_cosfilter<1>(nufftIn, nufftOut, 512, 2);
}


void nufft_cosfilter3(hls::stream< t_nufft_output_complex > nufftIn[3],
					 hls::stream< t_nufft_output_complex > nufftOut[3],
					 const int nL, const int m) {

#pragma HLS data_pack variable=nufftIn
#pragma HLS data_pack variable=nufftOut
#pragma HLS DATAFLOW
	nufft_cosfilter<3>(nufftIn, nufftOut, 512, 2);
}
