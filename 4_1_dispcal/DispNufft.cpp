
#include "DispNufft.h"
#include "hls_fpo.h"
#include "hls_math.h"
#include <math.h>
t_disp_scalar expo(t_disp_scalar dispVal){
	float a = hls::expf(dispVal.to_float());
	return t_disp_scalar(a);
}
void nufft_Disp (
		hls::stream<data_t>					&src_axi,
		hls::stream<t_disp_scalar>          &dispFilters){
#pragma HLS STREAM variable=dispFilters depth=7600 dim=1
#pragma HLS STREAM variable=src_axi depth=1520 dim=1
#pragma HLS INLINE
#pragma HLS DATAFLOW
#pragma HLS INTERFACE axis depth=7600 port=dispFilters
#pragma HLS INTERFACE axis depth=1520 port=src_axi

	hls::stream<t_disp_scalar>          dispFiltersInternal[5];//last for AA
//#pragma HLS RESOURCE variable=dispFiltersInternal core=RAM_2P_BRAM
#pragma HLS stream variable=dispFiltersInternal depth=1520

	const int climits[] = { 0, 512, 1024,1280,1408,1472,1504};
	const t_disp_scalar depPhsRatoio[] = {0.394429578492824,	0.557580824925606,	1.11513039602439,	2.23001066655451,	4.45801734656206,	8.89991391132190,	0.345847336083279};
	for(int i=0;i<1520;i++) {
#pragma HLS pipeline
	const t_disp_scalar  conV = 4.5;
	ap_fixed< 16,8> dA = 2.3;
	// Extracting DISP
	t_disp_scalar     depthScalarVal = t_disp_scalar(src_axi.read());
	t_disp_scalar dispVal;
	t_disp_scalar ratio;
	if(i>=climits[0]&&i<climits[1])ratio = depPhsRatoio[0];
	if(i>=climits[1]&&i<climits[2])ratio = depPhsRatoio[1];
	if(i>=climits[2]&&i<climits[3])ratio = depPhsRatoio[2];
	if(i>=climits[3]&&i<climits[4])ratio = depPhsRatoio[3];
	if(i>=climits[4]&&i<climits[5])ratio = depPhsRatoio[4];
	if(i>=climits[5]&&i<climits[6])ratio = depPhsRatoio[5];
	if(i>=climits[6])ratio = depPhsRatoio[6];

				dispVal = depthScalarVal / ratio;
				dispVal = dispVal * t_disp_scalar(0.3);
				dispVal = (dispVal * dispVal)>>1;
				dispVal = dispVal *-1;
				t_disp_scalar temp = dispVal;
				dispVal = expo(temp);
				//printf("%f %f %f\n",temp.to_float() ,a ,dispVal.to_float());

	for(int v=0;v<4;v++) {
	#pragma HLS loop_unroll
		t_disp_scalar displacement = ((conV - (v+1)) * dA - t_disp_scalar(0.5))*depthScalarVal;
					dispFiltersInternal[v+1].write(displacement);
	}
	dispFiltersInternal[0].write(dispVal);
	}
	for (int numR = 0; numR < 5; numR++) {
		for(int i=0;i<1520;i++) {
	#pragma HLS pipeline
		dispFilters.write(dispFiltersInternal[numR].read());
		//dispFilters[1].write(dispFiltersInternal[1].read());
		//dispFilters[2].write(dispFiltersInternal[2].read());
		//dispFilters[3].write(dispFiltersInternal[3].read());
	}
	}

}
