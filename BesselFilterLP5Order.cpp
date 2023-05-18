#include <math.h>
#include <stdint.h>
#include "BesselFilterLP5Order.h"


/**
 *  @brief Creates basic filter with cutoff frequency 1 rad/s and sample frequency 1 Hz
 */
BesselFilterLP5Order::BesselFilterLP5Order() : f_s{1.0}, cutoff_freq{2 * M_PI} {}

/**
 * @brief Constructs filter object
 * @param freq_sample sample frequency in Hz
 * @param cut_f cutoff filter frequency in Hz (is automatically transformed to rad/s)
 */
BesselFilterLP5Order::BesselFilterLP5Order(float freq_sample, float cut_f) : f_s{freq_sample} {
    set_cuttoff_frequency(cut_f);
}

BesselFilterLP5Order::~BesselFilterLP5Order() {}

/**
 * @brief Used to filtrate for one step in accordance to sample frequency
 * @param new_val new filter input value
 * @return filtered value y
 */
float BesselFilterLP5Order::step(float new_val) {
    float y = new_val * b[0] + x_prev[0] * b[1] + x_prev[1] * b[2] + 
            x_prev[2] * b[3] + x_prev[3] * b[4] + x_prev[4] * b[5] - 
            y_prev[0] * a[1] - y_prev[1] * a[2] - y_prev[2] * a[3] -
            y_prev[3] * a[4] - y_prev[4] * a[5];
    y /= a[0];

    /* delay for the x and y values */
    for (uint8_t i = (sizeof(x_prev) / sizeof(float)) - 1; i > 0; i--) {
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
void BesselFilterLP5Order::set_cuttoff_frequency(float new_cutoff) {
    /* changing to rad/s */
    cutoff_freq = new_cutoff * 2 * M_PI;
    /* recalculating of the coefficients */
    coefficients_calculating();
}

/**
 * @brief Method for getting the cutoff frequency value 
 * @return cutoff frequency in rad/s
 */
float BesselFilterLP5Order::get_cuttoff_frequency(void) {
    return cutoff_freq;
}

/**
 * @brief Method for setting the sample frequency value
 * Also method allows to recalculate automatically the filter coefficients 
 * @return cutoff frequency in Hz
 */
void BesselFilterLP5Order::set_sample_frequency(float new_f_s) {
    f_s = new_f_s;
    coefficients_calculating();
}

/**
 * @brief Getting the sample frequency value
 * @return sample frequency in Hz
 */
float BesselFilterLP5Order::get_sample_frequency(void) {
    return f_s;
}

/**
 * @brief Calculating the filter coefficients
 * Firstly, the analog filter coefficients according to cutoff frequency are obtained,
 * Finally, the digital filter coefficients are calculated
 * Coefficients equations were obtained with bilinear transformation of the analog filter transfer function
 */
void BesselFilterLP5Order::coefficients_calculating(void) {
    float tmp_a[5]; // temporary values to keep coefficients
    float a_0; // for divide operation
    float f_2, f_3, f_4, f_5; // for keeping the powers of sample frequency

    /*
    shifting the cuttoff frequency from 1 rad/s (lowpass mode), 
    caclulating the coefficients of the analog (continiuous) filter
    */
    for (uint8_t i = 0; i < (sizeof(tmp_a) / sizeof(float)); i++) {
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
