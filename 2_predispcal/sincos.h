#ifndef __HLS_SINCOS_FUNC_
#define __HLS_SINCOS_FUNC_

namespace hls {

  // ===================================================================================================================
  // Input and output interfaces
  template <int InputWidth> struct sincos_input {
    typedef cordic_inputs<InputWidth,CORDIC_F_SIN_COS,  hls::CORDIC_FORMAT_SIG_FRAC> in;
  };

  template <int OutputWidth> struct sincos_output {
    typedef cordic_outputs<OutputWidth,CORDIC_F_SIN_COS, hls::CORDIC_FORMAT_SIG_FRAC> out;
  };

  // ===================================================================================================================
  // SINCOS: Entry point function.
  // o Template parameters:
  //  - DataFormat       : Selects between unsigned fraction (with integer width of 1 bit) and unsigned integer formats
  //  - InputWidth       : Defines overall input data width
  //  - OutputWidth      : Defines overall output data width
  //  - RoundMode        : Selects the rounding mode to apply to the output data
  // o Arguments:
  //  - x                : Input data
  //  - sqrtX            : Square root of input data
  // o The internal CORDIC function applies its own rounding, therefore the interface
  //   ap_fixed types need not specify rounding and saturation modes
  template <
    int InputWidth,
    int OutputWidth,
    int RoundMode
    > void sincos(const typename sincos_input<InputWidth>::in &x,
                typename sincos_output<OutputWidth>::out &cordicOut) {

    Function_sincos_cordic:;

    cordic_base<
	CORDIC_F_SIN_COS,
      CORDIC_FALSE,
	  CORDIC_FORMAT_SIG_FRAC,
      CORDIC_FORMAT_RAD,
      InputWidth,
      OutputWidth,
      CORDIC_ITER_AUTO,
      CORDIC_PREC_AUTO,
      RoundMode,
      CORDIC_SCALE_NONE>(x, cordicOut);
  }

} // end namespace hls


#endif
