#pragma once

#include <numbers> // for pi
#include <concepts>

namespace BesselLowpassFilter {

/**
 *  *_cont variables for calculations of the filter's coefficients
    these variables are connected with analog lowpass Bessel filter with cutoff frequency 1 rad/s
 */

constexpr float b_0_cont = 1.0f;
constexpr float a_cont[] = {
  1.0f, 3.810701205349278f, 6.776673715676871f, 6.886367652423632f, 3.936283427035352f
};

template <typename T>
concept vals_to_filter = requires (T value) {
  value = value;
  value += value;
  value -= value;
  value + value;
  value - value;
  value / value;
  value * value;
  value < value;
  value > value;
  value <= value;
  value >= value;
};

template <typename F>
concept float_or_double = std::same_as<F, float> || std::same_as<F, double>;

/**
 * @brief Digital Bessel filter
 * 
 * Lowpass type, only 5th order
 */
template<vals_to_filter T, float_or_double F>
class BesselFilterLP5Order {
public:
  /**
   * @brief Used to filtrate for one step in accordance to sample frequency
   * @param new_val new filter input value
   * @return filtered value y
   */
  T step(T new_val) {
      T y = new_val * b[0] + x_prev[0] * b[1] + x_prev[1] * b[2] + 
        x_prev[2] * b[3] + x_prev[3] * b[4] + x_prev[4] * b[5] - 
        y_prev[0] * a[1] - y_prev[1] * a[2] - y_prev[2] * a[3] -
        y_prev[3] * a[4] - y_prev[4] * a[5];
      y /= a[0];

      /* delay for the x and y values */
      for (uint8_t i{(sizeof(x_prev) / sizeof(T)) - 1}; i > 0; i--) {
        x_prev[i] = x_prev[i - 1];
        y_prev[i] = y_prev[i - 1];
      }
      x_prev[0] = new_val;
      y_prev[0] = y;

      return y;
  }

  /**
   * @brief Setting new cutoff frequency with recalculation of the coefficients a and b
   * @param new_cutoff cutoff frequency in Hz
   */
  void set_cuttoff_frequency(F new_cutoff) {
    /* changing to rad/s */
    cutoff_freq = new_cutoff * 2 * std::numbers::pi;
    /* recalculating of the coefficients */
    coefficients_calculating();
  }

  /**
   * @brief Method for getting the cutoff frequency value 
   * @return cutoff frequency in rad/s
   */
  F get_cuttoff_frequency(void) {
    return cutoff_freq;
  }

  /**
   * @brief Setting new cutoff frequency with recalculation of the coefficients a and b
   * @param new_f_s         new sample rate frequency in Hz
   * @param new_cutoff_freq cutoff frequency in Hz
   */
  void set_cutoff_and_sample_frequencies(F new_f_s, F new_cutoff_freq) {
    f_s = new_f_s;
    cutoff_freq = new_cutoff_freq * 2 * std::numbers::pi;
    coefficients_calculating();
  }

  /**
   * @brief Method for setting the sample frequency value
   * Also method allows to recalculate automatically the filter coefficients 
   * @return cutoff frequency in Hz
   */
  void set_sample_frequency(F new_f_s) {
    f_s = new_f_s;
    coefficients_calculating();
  }

  /**
   * @brief Getting the sample frequency value
   * @return sample frequency in Hz
   */
  F get_sample_frequency(void) {
    return f_s;
  }

protected:
  /**
   * @brief Calculating the filter coefficients
   * Firstly, the analog filter coefficients according to cutoff frequency are obtained,
   * Finally, the digital filter coefficients are calculated
   * Coefficients equations were obtained with bilinear transformation of the analog filter transfer function
   */
  void coefficients_calculating(void) {
    F tmp_a[5]; // temporary values to keep coefficients
    F a_0; // for divide operation
    F f_2, f_3, f_4, f_5; // for keeping the powers of sample frequency

    /*
    shifting the cuttoff frequency from 1 rad/s (lowpass mode), 
    caclulating the coefficients of the analog (continiuous) filter
    */
    for (uint8_t i{0}; i < (sizeof(tmp_a) / sizeof(F)); i++) {
        tmp_a[i] = a_cont[i] / pow(cutoff_freq, 5 - i);
    }

    /* powers of the sample frequency */
    f_2 = pow(f_s, 2);
    f_3 = pow(f_s, 3);
    f_4 = pow(f_s, 4);
    f_5 = pow(f_s, 5);

    /* calculating the coefficients of the discrete filter (using the bilinear transformation) */
    a_0 = 1 + tmp_a[4] * 2 * f_s + tmp_a[3] * 4 * f_2 + tmp_a[2] * 8 * f_3 + tmp_a[1] * 16 * f_4 + tmp_a[0] * 32 * f_5;
    a[0] = 1.0;
    a[1] = (5 + 6 * tmp_a[4] * f_s + 4 * tmp_a[3] * f_2 - tmp_a[2] * 8 * f_3 - 48 * tmp_a[1] * f_4 - 160 * tmp_a[0] * f_5) / a_0;
    a[2] = (10 + 4 * tmp_a[4] * f_s - 8 * tmp_a[3] * f_2 - 16 * tmp_a[2] * f_3 + 32 * tmp_a[1] * f_4 + 320 * tmp_a[0] * f_5) / a_0;
    a[3] = (10 - 4 * tmp_a[4] * f_s - 8 * tmp_a[3] * f_2 + 16 * tmp_a[2] * f_3 + 32 * tmp_a[1] * f_4 - 320 * tmp_a[0] * f_5) / a_0;
    a[4] = (5 - 6 * tmp_a[4] * f_s + 4 * tmp_a[3] * f_2 + 8 * tmp_a[2] * f_3 - 48 * tmp_a[1] * f_4 + 160 * tmp_a[0] * f_5) / a_0;
    a[5] = (1 - 2 * tmp_a[4] * f_s + 4 * tmp_a[3] * f_2 - 8 * tmp_a[2] * f_3 + 16 * tmp_a[1] * f_4 - 32 * tmp_a[0] * f_5) / a_0; 

    b[0] = b_0_cont / a_0;
    b[1] = (b_0_cont * 5) / a_0;
    b[2] = (b_0_cont * 10) / a_0;
    b[3] = b[2];
    b[4] = b[1];
    b[5] = b[0];
  }

protected:
  F f_s{1.0}; // sample frequency in Hz
  F cutoff_freq{2 * std::numbers::pi}; // cutoff frequency in rad/s (default is 1 rad/s)
  F a[6], b[6]; // filter coefficients
  T y_prev[5] = {0,}, x_prev[5] = {0,}; // arrays for keeping previous x (input) and y (output) values
};

} /* namespace BesselLowpass */
