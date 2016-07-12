#include <math.h>
#include "common_type.h"

//pyrConstruct
//#include "../1_PyrConstruct/pyr.h"
////DispCheck
//#include "../2_DispCheck/cmpy_complex.h"
//#include "muxPyr_Disp.h"
////MedianFiler
//#include "../3_MedianFilter/median_filter.h"

void muxPyr_disp(
		hls::stream<t_pyr_complex> &pyrFilOut,
		hls::stream<t_input_complex> &sigs,
		hls::stream<t_input_complex> &sigRefs
);


void viewExpansionInterface(
		hls::stream<cmpxDataIn> &imgIn,
		hls::stream<t_disp_scalar> &prealign
){

#pragma HLS INTERFACE axis depth=2560 port=prealign
#pragma HLS INTERFACE axis depth=5120 port=imgIn
#pragma HLS DATA_PACK variable=prealign
#pragma HLS DATA_PACK variable=imgIn


	//PyrConstruct
	hls::stream<cmpxDataIn> imgin;
	hls::stream<t_pyr_complex> pyrFilOut;

#pragma HLS STREAM variable=imgin depth=512 dim=1
#pragma HLS DATA_PACK variable=pyrFilOut
#pragma HLS STREAM variable=pyrFilOut depth=15200 dim=1

	for (int i = 0; i < 10; i++) {
		for (int j = 0; j < 512; ++j) imgin.write(imgIn.read());
		pyrcon_top(imgin,pyrFilOut,512);


	}

#ifndef __SYNTHESIS__
	WritePlot("pyrOut.txt",pyrFilOut);
#endif

	hls::stream<t_input_complex> sigs;
	hls::stream<t_input_complex> sigRefs;

#pragma HLS STREAM variable=sigs depth=7600 dim=1
#pragma HLS DATA_PACK variable=sigs
#pragma HLS STREAM variable=sigRefs depth=7600 dim=1
#pragma HLS DATA_PACK variable=sigRefs
	muxPyr_disp(pyrFilOut,sigs,sigRefs);
	muxPyr_disp(pyrFilOut,sigs,sigRefs);
	muxPyr_disp(pyrFilOut,sigs,sigRefs);
	muxPyr_disp(pyrFilOut,sigs,sigRefs);
	muxPyr_disp(pyrFilOut,sigs,sigRefs);
	//Disp Check
	hls::stream<t_input_complex> sig;
	hls::stream<t_input_complex> sigRef;
	hls::stream<t_disp_scalar> Disp;
	hls::stream<t_output_complex> cmp;

	hls::stream<t_disp_scalar> PhaseDisp;
		int nL = 512;
		int nLExp = 1024;
		int width = 1024;
#pragma HLS STREAM variable=sig depth=1520 dim=1
#pragma HLS STREAM variable=sigRef depth=1520 dim=1
#pragma HLS STREAM variable=Disp depth=512 dim=1
#pragma HLS STREAM variable=cmp depth=7600 dim=1
#pragma HLS STREAM variable=PhaseDisp depth=7600 dim=1
#pragma HLS DATA_PACK variable=cmp
//		hls::stream<t_output_complex> cmpS;
//#pragma HLS STREAM variable=cmpS depth=1520 dim=1
//#pragma HLS DATA_PACK variable=cmpS
//		hls::stream<t_disp_scalar> Phaseout;
//#pragma HLS STREAM variable=Phaseout depth=2560 dim=1
//#pragma HLS DATA_PACK variable=Phaseout

	for (int i = 0; i < 5; ++i) {
#pragma HLS PIPELINE
		for (int j = 0; j < NLEN; ++j) {
			sig.write(sigs.read());
			sigRef.write(sigRefs.read());
		}
		for (int j = 0; j < 512; ++j) {
			Disp.write(prealign.read());
				}
		cmpy_complex_top(sig,sigRef,Disp,cmp,PhaseDisp,nL,nLExp,nL,-1.0*nL/width);
}
#ifndef __SYNTHESIS__
	WritePlot("DispOutCmpy.txt",cmp);
	WritePlotsca("DispOutDepth.txt",PhaseDisp);
#endif

	hls::stream<data_t> PreDisp;
#pragma HLS DATA_PACK variable=PreDisp

	hls::stream<data_t> PostDisp;
#pragma HLS DATA_PACK variable=PostDisp
#pragma HLS STREAM variable=PreDisp depth=7600 dim=1
#pragma HLS STREAM variable=PostDisp depth=7600 dim=1


for (int i = 0;  i < 7600;  i++) {
	data_t temp = PhaseDisp.read();
	PreDisp.write(temp);
}
	median_strm(1520,5,PreDisp,PostDisp);
	//WritePlotsca("PostDepth.txt",PostDisp);


#ifndef __SYNTHESIS__
	FILE * fo = fopen("mdf_out.txt", " wb");
	for (int j = 0; j <5; ++j) {

	for( int i=0;i<1520;i++)
		{
			data_t val2 = PostDisp.read();
			fprintf(fo, "%.8f  , ", val2.to_double());
				}
	fprintf(fo," \n");
	}

#endif
}


