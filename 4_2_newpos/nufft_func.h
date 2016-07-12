#ifndef __NUFFT_FUNC__
#define __NUFFT_FUNC__

#include "coeff.h"

#define   Q    4
#define   M    2

//      const _Tp __r = _M_real * __z.real() - _M_imag * __z.imag();
//      _M_imag = _M_real * __z.imag() + _M_imag * __z.real();
//      _M_real = __r;
//      return *this;
template <typename T>
std::complex<T> cmpMulPipeline(const std::complex<T> &A, const std::complex<T> &B) {
	T real = A.real() * B.real() -  A.imag() * B.imag();
	T imag = A.real() * B.imag() +  A.imag() * B.real();

#pragma HLS RESOURCE variable=real core=Mul2S
#pragma HLS RESOURCE variable=imag core=Mul2S
	return std::complex<T> (real, imag);
}


template<int C, int MAXWIDTH>
void nufft_top(hls::stream<t_input_complex>  &sig,
			   hls::stream<t_disp_scalar>    &dispFilter,
			   hls::stream<t_nufft_output_complex> &sigStreamOut,
			   const int nL,
			   int mask) {

	//const int nL = 128;
	//const int Kx = 0;
#pragma HLS DATAFLOW

#pragma HLS data_pack variable=c0
#pragma HLS data_pack variable=c1
#pragma HLS data_pack variable=c2
#pragma HLS data_pack variable=c3

#pragma HLS data_pack variable=sig
#pragma HLS data_pack variable=sigStreamOut

//#define DEBUG

	const t_output_scalar pi2InvNlM = M_PI * 2 / (M * nL);
	const t_output_scalar piInvM = M_PI / M;

	t_nufft_output_complex sigMod0[MAXWIDTH * M /4];
	t_nufft_output_complex sigMod1[MAXWIDTH * M /4];
	t_nufft_output_complex sigMod2[MAXWIDTH * M /4];
	t_nufft_output_complex sigMod3[MAXWIDTH * M /4];

//#pragma HLS ARRAY_PARTITION variable=sigMod0 block factor=2 dim=1
#pragma HLS data_pack       variable=sigMod0
//#pragma HLS RESOURCE        variable=sigMod0 core=RAM_2P

//#pragma HLS ARRAY_PARTITION variable=sigMod1 block factor=2 dim=1
#pragma HLS data_pack       variable=sigMod1
//#pragma HLS RESOURCE        variable=sigMod1 core=RAM_2P

//#pragma HLS ARRAY_PARTITION variable=sigMod2 block factor=2 dim=1
#pragma HLS data_pack       variable=sigMod2
//#pragma HLS RESOURCE        variable=sigMod2 core=RAM_2P

//#pragma HLS ARRAY_PARTITION variable=sigMod3 block factor=2 dim=1
#pragma HLS data_pack       variable=sigMod3
//#pragma HLS RESOURCE        variable=sigMod3 core=RAM_2P


	for(int i=0;i<nL * M/4;i++ ) {
#pragma HLS loop_tripcount min=16 avg=128 max=512
#pragma HLS PIPELINE
		const t_nufft_output_complex  zero(0,0);
		//zero.real() = 0;
		//zero.imag() = 0;

			sigMod0[i] = zero;
			sigMod1[i] = zero;
			sigMod2[i] = zero;
			sigMod3[i] = zero;

	}
	WAVELETLOOP: for(int i=0;i<nL;i++ ) {
#pragma HLS INLINE
#pragma HLS loop_tripcount min=16 avg=128 max=512
//#pragma HLS PIPELINE

		t_disp_scalar disp = dispFilter.read();
		t_disp_scalar newIndex =  (disp + i + 1) * M  - ((t_disp_scalar)Q-1)/2;
		if (newIndex <0) newIndex = 0;
	//	if (newIndex >= (nL*M)-4) newIndex = nL*M - 4;

		t_disp_scalar newIndex_floor = (int)newIndex;

		ap_int<12> frac_index = ((newIndex - newIndex_floor) << 10) + t_disp_scalar(0.5);

		t_nufft_output_complex X0 = c0[frac_index];
		t_nufft_output_complex X1 = c1[frac_index];
		t_nufft_output_complex X2 = c2[frac_index];
		t_nufft_output_complex X3 = c3[frac_index];


//#pragma HLS unroll
			int nIdx = newIndex_floor;

			t_nufft_output_complex pyrVal = sig.read();

//			if (c==0 && nL==256)
//				printf("pyrval: %.8f %.8f\n", pyrVal.real().to_float(), pyrVal.imag().to_float());

			t_nufft_output_complex pyrValMulX0 = cmpMulPipeline< t_nufft_output_scalar>(pyrVal , X0);
			t_nufft_output_complex pyrValMulX1 = cmpMulPipeline< t_nufft_output_scalar>(pyrVal , X1);
			t_nufft_output_complex pyrValMulX2 = cmpMulPipeline< t_nufft_output_scalar>(pyrVal , X2);
			t_nufft_output_complex pyrValMulX3 = cmpMulPipeline< t_nufft_output_scalar>(pyrVal , X3);


#pragma HLS data_pack variable=pyrVal


#pragma HLS data_pack variable=pyrValMulX0
#pragma HLS RESOURCE variable=pyrValMulX0 core=Mul3S
#pragma HLS data_pack variable=pyrValMulX1
#pragma HLS RESOURCE  variable=pyrValMulX1 core=Mul3S
#pragma HLS data_pack variable=pyrValMulX2
#pragma HLS RESOURCE  variable=pyrValMulX2 core=Mul3S
#pragma HLS data_pack variable=pyrValMulX3
#pragma HLS RESOURCE  variable=pyrValMulX3 core=Mul3S
//#pragma HLS RESOURCE  variable=pyrValMulX3 core=Mul3S
			unsigned short K = nIdx >> 2;
			unsigned short Kleft = nIdx & 3;

			unsigned short nIdx0 = K + (Kleft>0);
			unsigned short nIdx1 = K + (Kleft>1);
			unsigned short nIdx2 = K + (Kleft>2);
			unsigned short nIdx3 = K;

			nIdx0 &= mask;
			nIdx1 &= mask;
			nIdx2 &= mask;
			nIdx3 &= mask;
			t_output_complex pyrValX[4];
#pragma HLS ARRAY_PARTITION variable=pyrValX complete dim=1
#pragma HLS data_pack variable=pyrValX

			if (Kleft == 0) {
				pyrValX[0] = pyrValMulX0;
				pyrValX[1] = pyrValMulX1;
				pyrValX[2] = pyrValMulX2;
				pyrValX[3] = pyrValMulX3;
			}
			if (Kleft == 1) {
				pyrValX[0] = pyrValMulX3;
				pyrValX[1] = pyrValMulX0;
				pyrValX[2] = pyrValMulX1;
				pyrValX[3] = pyrValMulX2;
			}
			if (Kleft == 2) {
				pyrValX[0] = pyrValMulX2;
				pyrValX[1] = pyrValMulX3;
				pyrValX[2] = pyrValMulX0;
				pyrValX[3] = pyrValMulX1;
			}
			if (Kleft == 3) {
				pyrValX[0] = pyrValMulX1;
				pyrValX[1] = pyrValMulX2;
				pyrValX[2] = pyrValMulX3;
				pyrValX[3] = pyrValMulX0;
			}




			sigMod0[ nIdx0 ] += pyrValX[0];
			sigMod1[ nIdx1 ] += pyrValX[1];
			sigMod2[ nIdx2 ] += pyrValX[2];
			sigMod3[ nIdx3 ] += pyrValX[3];

	}
	for(int i=0;i<nL * M/4;i++ ) {
#pragma HLS loop_tripcount min=16 avg=128 max=1024

#pragma HLS PIPELINE II=4

//#pragma HLS unroll
			sigStreamOut.write( sigMod0[i]);
			sigStreamOut.write( sigMod1[i]);
			sigStreamOut.write( sigMod2[i]);
			sigStreamOut.write( sigMod3[i]);

	}
}


template<int C, int MAXWIDTH>
void nufft_top(hls::stream<t_input_complex>  &sig,
			   hls::stream<t_disp_scalar>    &dispFilter,
			   hls::stream<t_nufft_output_complex> &sigStreamOut,
			   const int nL,
			   const int Kx,
			   int mask) {

	//const int nL = 128;
	//const int Kx = 0;
//#pragma HLS DATAFLOW

#pragma HLS data_pack variable=cL0
#pragma HLS data_pack variable=cL1
#pragma HLS data_pack variable=cL2
#pragma HLS data_pack variable=cL3

#pragma HLS data_pack variable=sig
#pragma HLS data_pack variable=sigStreamOut

//#define DEBUG

	const t_output_scalar pi2InvNlM = M_PI * 2 / (M * nL);
	const t_output_scalar piInvM = M_PI / M;

	t_nufft_output_complex sigMod0[MAXWIDTH * M /4];
	t_nufft_output_complex sigMod1[MAXWIDTH * M /4];
	t_nufft_output_complex sigMod2[MAXWIDTH * M /4];
	t_nufft_output_complex sigMod3[MAXWIDTH * M /4];

#pragma HLS ARRAY_PARTITION variable=sigMod0 block factor=3 dim=1
#pragma HLS data_pack       variable=sigMod0
#pragma HLS RESOURCE        variable=sigMod0 core=RAM_2P

#pragma HLS ARRAY_PARTITION variable=sigMod1 block factor=3 dim=1
#pragma HLS data_pack       variable=sigMod1
#pragma HLS RESOURCE        variable=sigMod1 core=RAM_2P

	#pragma HLS ARRAY_PARTITION variable=sigMod2 block factor=3 dim=1
#pragma HLS data_pack       variable=sigMod2
#pragma HLS RESOURCE        variable=sigMod2 core=RAM_2P

#pragma HLS ARRAY_PARTITION variable=sigMod3 block factor=3 dim=1
#pragma HLS data_pack       variable=sigMod3
#pragma HLS RESOURCE        variable=sigMod3 core=RAM_2P


	for(int i=0;i<nL * M/4;i++ ) {
#pragma HLS loop_tripcount min=16 avg=128 max=512
#pragma HLS PIPELINE
		const t_nufft_output_complex  zero(0,0);
		//zero.real() = 0;
		//zero.imag() = 0;

			sigMod0[i] = zero;
			sigMod1[i] = zero;
			sigMod2[i] = zero;
			sigMod3[i] = zero;

	}
	WAVELETLOOP: for(int i=0;i<nL;i++ ) {
#pragma HLS INLINE
#pragma HLS loop_tripcount min=16 avg=128 max=512
#pragma HLS PIPELINE II=4

		t_disp_scalar disp = dispFilter.read();
		t_disp_scalar newIndex =  (disp + i + 1) * M  - ((t_disp_scalar)Q-1)/2;
		if (newIndex <0) newIndex = 0;
	//	if (newIndex >= (nL*M)-4) newIndex = nL*M - 4;

		t_disp_scalar newIndex_floor = (int)newIndex;

		ap_int<12> frac_index = ((newIndex - newIndex_floor) << 10) + t_disp_scalar(0.5);

		t_nufft_output_complex X0 = cL0[Kx][frac_index];
		t_nufft_output_complex X1 = cL1[Kx][frac_index];
		t_nufft_output_complex X2 = cL2[Kx][frac_index];
		t_nufft_output_complex X3 = cL3[Kx][frac_index];


//#pragma HLS unroll
			int nIdx = newIndex_floor;

			t_nufft_output_complex pyrVal = sig.read();

//			if (c==0 && nL==256)
//				printf("pyrval: %.8f %.8f\n", pyrVal.real().to_float(), pyrVal.imag().to_float());

			t_nufft_output_complex pyrValMulX0 = pyrVal * X0;
			t_nufft_output_complex pyrValMulX1 = pyrVal * X1;
			t_nufft_output_complex pyrValMulX2 = pyrVal * X2;
			t_nufft_output_complex pyrValMulX3 = pyrVal * X3;


#pragma HLS data_pack variable=pyrVal


#pragma HLS data_pack variable=pyrValMulX0
#pragma HLS RESOURCE variable=pyrValMulX0 core=Mul3S
#pragma HLS data_pack variable=pyrValMulX1
#pragma HLS RESOURCE  variable=pyrValMulX1 core=Mul3S
#pragma HLS data_pack variable=pyrValMulX2
#pragma HLS RESOURCE  variable=pyrValMulX2 core=Mul3S
#pragma HLS data_pack variable=pyrValMulX3
#pragma HLS RESOURCE  variable=pyrValMulX3 core=Mul3S

			unsigned short K = nIdx >> 2;
			unsigned short Kleft = nIdx & 3;

			unsigned short nIdx0 = K + (Kleft>0);
			unsigned short nIdx1 = K + (Kleft>1);
			unsigned short nIdx2 = K + (Kleft>2);
			unsigned short nIdx3 = K;

			nIdx0 &= mask;
			nIdx1 &= mask;
			nIdx2 &= mask;
			nIdx3 &= mask;
			t_output_complex pyrValX[4];
#pragma HLS data_pack variable=pyrValX
#pragma HLS ARRAY_PARTITION variable=pyrValX complete dim=1

//			pyrValX[0] = Kleft == 0 ?
			if (Kleft == 0) {
				pyrValX[0] = pyrValMulX0;
				pyrValX[1] = pyrValMulX1;
				pyrValX[2] = pyrValMulX2;
				pyrValX[3] = pyrValMulX3;
			}
			if (Kleft == 1) {
				pyrValX[0] = pyrValMulX3;
				pyrValX[1] = pyrValMulX0;
				pyrValX[2] = pyrValMulX1;
				pyrValX[3] = pyrValMulX2;
			}
			if (Kleft == 2) {
				pyrValX[0] = pyrValMulX2;
				pyrValX[1] = pyrValMulX3;
				pyrValX[2] = pyrValMulX0;
				pyrValX[3] = pyrValMulX1;
			}
			if (Kleft == 3) {
				pyrValX[0] = pyrValMulX1;
				pyrValX[1] = pyrValMulX2;
				pyrValX[2] = pyrValMulX3;
				pyrValX[3] = pyrValMulX0;
			}

			sigMod0[ nIdx0 ] += pyrValX[0];
			sigMod1[ nIdx1 ] += pyrValX[1];
			sigMod2[ nIdx2 ] += pyrValX[2];
			sigMod3[ nIdx3 ] += pyrValX[3];



	}
	for(int i=0;i<nL * M/4;i++ ) {
#pragma HLS loop_tripcount min=16 avg=128 max=1024

#pragma HLS PIPELINE II=4

			sigStreamOut.write( sigMod0[i]);
			sigStreamOut.write( sigMod1[i]);
			sigStreamOut.write( sigMod2[i]);
			sigStreamOut.write( sigMod3[i]);

	}
}

#endif
