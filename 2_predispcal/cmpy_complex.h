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
#ifndef CMPY_COMPLEX_H
#define CMPY_COMPLEX_H

#include "hls_dsp.h"
//#include "../TestFlow/common_type.h"

//#define NLEN 1024
#define NLEN 1520
//Configuration parameters for this instance.
typedef hls::CmpyThreeMult ARCH;
//ROUND_MODE Values : AP_RND, AP_RND_ZERO, AP_RND_MIN_INF, AP_RND_INF, AP_RND_CONV, AP_TRN, AP_TRN_ZERO
//OVERFLOW_MODE Values : AP_SAT, AP_SAT_ZERO, AP_SAT_SYM, AP_WRAP, AP_WRAP_SM

const int INPUT_WIDTH                = 16;//10;
const int INPUT_INTEGER_BITS         = 2;//2;
const int INPUT_ROUND_MODE           = 5;//AP_TRN, etc;
const int INPUT_OVERFLOW_MODE        = 3;//AP_WRAP, etc;
const int INPUT_SATURATION_BITS      = 0;

const int OUTPUT_WIDTH               = 20;
const int OUTPUT_INTEGER_BITS        = 3;//2
const int OUTPUT_ROUND_MODE          = 0;//AP_TRN, etc;
const int OUTPUT_OVERFLOW_MODE       = 3;//AP_WRAP, etc;
const int OUTPUT_SATURATION_BITS     = 0;

//Interface types
typedef ap_fixed<
  INPUT_WIDTH,
  INPUT_INTEGER_BITS,
  (ap_q_mode)INPUT_ROUND_MODE,
  (ap_o_mode)INPUT_OVERFLOW_MODE,
  INPUT_SATURATION_BITS> t_input_scalar;


const int DISP_WIDTH                = 24;
const int DISP_INTEGER_BITS         = 12;
const int DISP_ROUND_MODE           = 5;//AP_TRN, etc;
const int DISP_OVERFLOW_MODE        = 3;//AP_WRAP, etc;
const int DISP_SATURATION_BITS      = 0;

//Interface types
typedef ap_fixed<
		DISP_WIDTH,
		DISP_INTEGER_BITS,
  (ap_q_mode)DISP_ROUND_MODE,
  (ap_o_mode)DISP_OVERFLOW_MODE,
  DISP_SATURATION_BITS> t_disp_scalar;


typedef ap_fixed<
  OUTPUT_WIDTH,
  OUTPUT_INTEGER_BITS,
  (ap_q_mode)OUTPUT_ROUND_MODE,
  (ap_o_mode)OUTPUT_OVERFLOW_MODE,
  OUTPUT_SATURATION_BITS> t_output_scalar;

typedef std::complex< t_input_scalar > t_input_complex;
typedef std::complex< t_output_scalar > t_output_complex;

struct input{
	t_input_complex data;

	//bool last;
};
struct output{
	t_output_complex data;
	//bool last;
};
void cmpy_complex_top(
		hls::stream<t_input_complex> &sig,
		hls::stream<t_input_complex> &sigRef,
		hls::stream<t_disp_scalar> &prealign,
		hls::stream<t_output_complex> &cmp,
		hls::stream<t_disp_scalar> &PhaseDisp
   );

#endif
