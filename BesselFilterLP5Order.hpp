#pragma once

/**
 *  *_cont variables for calculations of the filter's coefficients
    these variables are connected with analog lowpass Bessel filter with cutoff frequency 1 rad/s
 */
constexpr float b_0_cont = 1.0; 
constexpr float a_cont[] = {1.0, 3.810701205349278, 6.776673715676871, 6.886367652423632, 3.936283427035352};

/**
 * @brief Class of the digital Bessel filter
    Lowpass type, only 5th order
 */
class BesselFilterLP5Order {
    public:
        BesselFilterLP5Order();
        BesselFilterLP5Order(float freq_sample, float cut_f);
        ~BesselFilterLP5Order();
        float step(float new_val);
        void set_cuttoff_frequency(float new_cutoff);
        float get_cuttoff_frequency(void);
        void set_sample_frequency(float new_f_s);
        float get_sample_frequency(void);
    protected:
        void coefficients_calculating(void);
        float f_s, cutoff_freq; // sample frequency in Hz and cutoff frequency in rad/s 
        float a[6], b[6]; // filter coefficients
        float y_prev[5], x_prev[5]; // arrays for keeping previous x (input) and y (output) values
};
