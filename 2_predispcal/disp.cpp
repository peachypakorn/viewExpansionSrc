/*****************************************************************************
 *
 *     Author: Xilinx, Inc.
 *
 *     This text contains proprietary, confidential information of
 *     Xilinx, Inc. , is distributed by under license from Xilinx,
 *     Inc., and may be used, copied and/or disclosed only pursuant to
 *     the terms of a valid license agreement with Xilinx, Inc.
 *
 *     XILINX IS PROVIDING THIS DESIGN, CODE, OR INFORMATION "AS IS"
 *     AS A COURTESY TO YOU, SOLELY FOR USE IN DEVELOPING PROGRAMS AND
 *     SOLUTIONS FOR XILINX DEVICES.  BY PROVIDING THIS DESIGN, CODE,
 *     OR INFORMATION AS ONE POSSIBLE IMPLEMENTATION OF THIS FEATURE,
 *     APPLICATION OR STANDARD, XILINX IS MAKING NO REPRESENTATION
 *     THAT THIS IMPLEMENTATION IS FREE FROM ANY CLAIMS OF INFRINGEMENT,
 *     AND YOU ARE RESPONSIBLE FOR OBTAINING ANY RIGHTS YOU MAY REQUIRE
 *     FOR YOUR IMPLEMENTATION.  XILINX EXPRESSLY DISCLAIMS ANY
 *     WARRANTY WHATSOEVER WITH RESPECT TO THE ADEQUACY OF THE
 *     IMPLEMENTATION, INCLUDING BUT NOT LIMITED TO ANY WARRANTIES OR
 *     REPRESENTATIONS THAT THIS IMPLEMENTATION IS FREE FROM CLAIMS OF
 *     INFRINGEMENT, IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *     FOR A PARTICULAR PURPOSE.
 *
 *     Xilinx products are not intended for use in life support appliances,
 *     devices, or systems. Use in such applications is expressly prohibited.
 *
 *     (c) Copyright 2014 Xilinx Inc.
 *     All rights reserved.
 *
 *****************************************************************************/

#include "cmpy_complex.h"
#include "hls_dsp.h"
#include "sincos.h"
#include "fxp_sqrt.h"

const int PhaseFormat = hls::CORDIC_FORMAT_RAD;
const int InputWidth = 16;
const int OutputWidth = 16;
const t_disp_scalar depPhsRatoio[] = {0.394429578492824,	0.557580824925606,	1.11513039602439,	2.23001066655451,	4.45801734656206,	8.89991391132190,	0.345847336083279};

//const int RoundMode = hls::CORDIC_ROUND_TRUNCATE;
const int RoundMode = hls::CORDIC_ROUND_POS_NEG_INF;
//  Complex *sig, Complex *sigRef, float *prealign, float *cmpR, float *cmpI,
//        int nL, int height, int nLExp, float factor, bool initialized, cudaStream_t stream



// The top-level function to synthesize

void atan2_top(const hls::atan2_input<InputWidth>::cartesian &x,
               hls::atan2_output<OutputWidth>::phase &atanX){

  // Call arctan function
  hls::atan2<PhaseFormat,InputWidth,OutputWidth,RoundMode>(x, atanX);
}

const int DataFormat = hls::CORDIC_FORMAT_USIG_FRAC;
//void sqrt_top(const hls::sqrt_input<InputWidth, DataFormat>::in &x,
//               hls::sqrt_output<OutputWidth, DataFormat>::out &sqrtx){


// The top-level function to synthesize
//
void sqrt_top2(const hls::sqrt_input<InputWidth, DataFormat>::in &x,
              hls::sqrt_output<OutputWidth, DataFormat>::out &sqrtX){
  // Call square root function
  hls::sqrt<DataFormat,InputWidth,OutputWidth,RoundMode>(x, sqrtX);
}

void sqrt_top( const t_output_scalar &x,
		             t_output_scalar  &xout) {
#pragma HLS PIPELINE
	const ap_ufixed< t_output_scalar::width, t_output_scalar::iwidth> xu = x;
	ap_ufixed< t_output_scalar::width, t_output_scalar::iwidth> xuout ;


//	hls::sqrt_input<InputWidth, DataFormat>::in  x_unsigned;
//	hls::sqrt_output<OutputWidth, DataFormat>::out  xout_unsigned;
	//x_unsigned.in = x;
	//sqrt_top2( x_unsigned, xout_unsigned);
	//xout = xout_unsigned.out;
	fxp_sqrt(xuout, xu);
	xout = xuout;
}
typedef ap_fixed<16,2,(ap_q_mode)6,(ap_o_mode)3,0> in;
ap_fixed<16,2> z8 = 0.0008;
ap_fixed<16,2> mz8 = -0.0008;
ap_fixed<16,2> z199 = 1.99;
ap_fixed<16,2> mz199 = -1.99;
template <class A, class B>
void myatan2( const A &xin,
		            B &yout

					)
{
	 hls::atan2_input<InputWidth>::cartesian  x;
	 in tmpR = xin.real();
	 in tmpI = xin.imag();

	 for (int i = 0; i < 4; ++i) {
		 	 if(tmpR>z199 || tmpR <mz199 || tmpI>z199 || tmpI <mz199){
		 		 tmpR = tmpR/2;
		 		 tmpI = tmpI/2;
		 	 }
		 	 else if((tmpR>0&&tmpR<z8)||(tmpR<0&&tmpR>mz8)||(tmpI>0&&tmpI<z8)||(tmpI<0&&tmpI>mz8)){
		 		tmpR = tmpR<<3;
		 		tmpI = tmpI<<3;
		 	 }
		}
//	 if(xin.real()<zeroP8 &&xin.imag()<zeroP8&&xin.real()>-zeroP8 &&xin.imag()>-zeroP8){
//		 x.cartesian.real() = in(xin.real())*2;
//		 x.cartesian.imag() = in(xin.imag())*2;
//	 }
	 //else{
		 x.cartesian.real() = tmpR;//in(xin.real());
		  x.cartesian.imag() =tmpI;//in(xin.imag()) ;

	 //}
		  //x.  sigRef[i];
		  //refAtans[i]

     hls::atan2_output<OutputWidth>::phase phase;
	 atan2_top(x, phase);
//     printf("%f %f %f\n ", x.cartesian.real().to_float(),
//    		 	 	 	   x.cartesian.imag().to_float(),
//				           phase.phase.to_float());
	 yout = phase.phase;
}


void sincos(const t_output_scalar &t, t_output_complex &scout ) {
#pragma HLS PIPELINE
	hls::sincos_input<InputWidth>::in  x;
	hls::sincos_output<OutputWidth>::out  xout;

	x.phase = t;


	//	hls::cordic_base()
	hls::sincos<InputWidth,OutputWidth,RoundMode>( x, xout);
	scout = xout.cartesian;

}

template <class A>
A  myreaminder(A x, A y) {
#pragma HLS PIPELINE
	A xdy = x/y;
	const A roundC = 0.5f;

	return x - int(xdy + roundC) * y;

}

const int Kset = 7;
const int limits[] = { 512, 512,256,128,64,32,16};
const int climits[] = { 0, 512, 1024,1280,1408,1472,1504};
const int llimits[] = { 9,9,8,7,6,5,4};
#define DISP 512


void cmpy_complex_top(
		hls::stream<t_input_complex> &sig,
		hls::stream<t_input_complex> &sigRef,
		hls::stream<t_disp_scalar> &prealign,
		hls::stream<t_output_complex> &cmp,
		hls::stream<t_disp_scalar> &PhaseDisp
             ){

#pragma HLS DATA_PACK variable=cmp
#pragma HLS DATA_PACK variable=sig
#pragma HLS data_pack variable=sigRef
#pragma HLS data_pack variable=prealign

#pragma HLS DATAFLOW
#pragma HLS INTERFACE axis port=sig
#pragma HLS INTERFACE axis port=sigRef
#pragma HLS INTERFACE axis port=prealign
#pragma HLS INTERFACE axis port=cmp
#pragma HLS INTERFACE axis port=PhaseDisp

#pragma HLS STREAM variable=sig depth=1520 dim=1
#pragma HLS STREAM variable=sigRef depth=1520 dim=1
#pragma HLS STREAM variable=prealign depth=1520 dim=1
#pragma HLS STREAM variable=cmp	depth=1520 dim=1
#pragma HLS STREAM variable=PhaseDisp depth=1520 dim=1

	t_disp_scalar disp[1520];
#pragma HLS STREAM variable=disp depth=1520 dim=1

	t_output_scalar refAtans[NLEN];
#pragma HLS RESOURCE variable=refAtans core=RAM_2P_BRAM

	for(int i=0;i<1520;i++) {
#pragma HLS PIPELINE
		disp[i] = prealign.read();
	}
	//start to find arctan of RefPic
	for(int i=0;i<NLEN;i++) {
#pragma HLS PIPELINE
	  t_input_complex temp = sigRef.read();
	 if(i>1503)refAtans[i] =0;
	 else myatan2(temp, refAtans[i]);//,temp.real().to_float(),temp.imag().to_float() );

  }
#ifndef __SYNTHESIS__
	{
		FILE * fo = fopen("PharseRef_out_1.txt", "wb");
		for(int i=0;i<1520;i++) {
			fprintf(fo, "%.8f \n",refAtans[i].to_float() );
		}
		fclose(fo);
	}
#endif

	//for debugging
	//t_output_scalar CheckAngle[NLEN];
	//t_output_complex CheckSinCos[NLEN];
	t_output_scalar CheckComplex[NLEN];

	int l = 0;
	int i = 0;
  for(int x=0;x<NLEN;x++ ){
#pragma HLS PIPELINE
	  int nlimit = limits[l];
	  int half =nlimit>>1;
	  int idx;
	  int limit = limits[l];
	  int cidx  = climits[l];
	  t_output_scalar angle;

	  if(i<half) idx = i;
	  else idx = 512 - nlimit + i;
	  //askkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk

	  t_input_scalar factor2 =  ap_ufixed<32,10>(limit)>>9;
	  t_disp_scalar pa = disp[x] * factor2 ;
	  t_disp_scalar xRef = (x + pa);
	  pa = xRef - x;
	 // printf("%d  %f\n",x , xRef.to_float() );
	 	  //run pixel by pixel so read from x
	  t_input_complex s = sig.read();
	  myatan2(s, angle);//,s.real().to_float(),s.imag().to_float());
	  if(x>1503){
		  angle = 0;
		//  xRef = 1;
	  }

	  t_disp_scalar dRes = refAtans[ xRef] - angle;
	  //printf("%d  %f  %f  %f  %f %f %f\n",x ,dRes.to_float(), refAtans[ xRef].to_float(),angle.to_float(),xRef.to_float(),factor2.to_float(),disp[x].to_float());
	  //dRes = 2*remainderf(dRes / M_PI + 2, 2);  INCOMPLETED
	  const t_output_scalar  mypi = M_PI;
	  //dRes =   ( (dRes / mypi ) ) ;  // To be implement, INCOMPLETED
	  const t_disp_scalar  val2 = 2;

	  //t_output_scalar tmp = dRes/mypi + val2;
	  t_disp_scalar tmp = dRes/mypi + val2;
	  tmp = val2 * myreaminder(tmp, val2);

	 printf("%f %f\n", tmp.to_float(), dRes.to_float());

	  // sincospif(dRes - pa, &aIm, &aRe);
	  t_output_complex sincosOut;
	 //PhaseDisp[x] =tmp;
	 tmp = (tmp*depPhsRatoio[l]) - disp[x];
	 if(tmp>mypi||tmp<(mypi*-1)){
		 t_output_scalar ovf =  tmp/mypi;
		 int mphs = (int)(ovf);
	 			  dRes = (ovf - mphs)*mypi;
	 		  }
	 	 else dRes = tmp;

	// printf("%f  %f \n" , tmp.to_float() , dRes.to_float());
	 // tmp = tmp - t_disp_scalar(pa/factor2);
	  ///////////////////////////////////////dRes = dRes - pa;
	 // a = dRes;

	  sincos( dRes, sincosOut);
	  //CheckSinCos[x] = sincosOut;
//	  printf("%f %f %f\n", dRes.to_float(),
//			  sincosOut.real().to_float(), sincosOut.imag().to_float());
	  // compute len

bool check = false;
const t_output_scalar  zpn = 0.1;
const t_output_scalar  mzpn = -0.1;

	  t_output_scalar tr ;
	  t_output_scalar ti;
  if((s.real()<zpn&&s.real()>mzpn)&&(s.imag()<zpn&&s.imag()>mzpn)) check = true;

  if(check)
	  {
		  tr = s.real()<<3;
		  ti = s.imag()<<3;
	  }
	  else {
		  tr = s.real();
		  ti = s.imag();
	  }
	  t_output_scalar val = tr*tr + ti *ti;
//	 if(x>1504){
//		 //val = (s.real()*s.real()+s.imag()*s.imag());
//		 val = val << 1;
//	 }

	  t_output_scalar lenT;
	  sqrt_top(val, lenT);
	  t_output_scalar len = lenT;
	 // if(x>1504) len = len <<2;
	  if(check) len = len >>3;


	 // printf("%f %f\n", val.to_float(), len.to_float());
	  t_output_complex sincosOut2;
	  sincosOut2.real() = s.real();//sincosOut.real() * len;
	  sincosOut2.imag() = s.imag();//sincosOut.imag() * len;
	  //otmp.data = sincosOut;
	  if(x>1503){
			  sincosOut2.real() = s.real();
			  sincosOut2.imag() = s.imag();
		  }
	  cmp.write(sincosOut2);
	  t_disp_scalar tmp2 = tmp;
	  PhaseDisp.write(tmp2);
	  i++;
	  l += (nlimit == i);
	  i = (nlimit == i)?0:i;
  }
#if 0
	{
		FILE * fo = fopen("Pharse_out_Angle.txt", "wb");
		for(int i=0;i<512;i++) {
			fprintf(fo, "%.8f \n",CheckAngle[i].to_float() );
		}
		fclose(fo);
	}
	{
			FILE * fo = fopen("SincosCheck.txt", "wb");
			for(int i=0;i<512;i++) {
				fprintf(fo, "%.8f %.8f  \n",CheckSinCos[i].real().to_float(),CheckSinCos[i].imag().to_float() );
			}
			fclose(fo);
		}
#endif
#ifndef __SYNTHESIS__
	{
			FILE * fo = fopen("CheckLength.txt", "wb");
			for(int i=0;i<1520;i++) {
				fprintf(fo, "%.8f \n",CheckComplex[i].to_float() );
			}
			fclose(fo);
		}
#endif


} // end of function cmpy_complex_top


