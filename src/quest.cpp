#include "quest.h"
using namespace godot;


/*
Note before everything, this implementation is directly tanslation from the python code 
No optimization has done. I did this way for two reason, although first this algorithm should 
not cause a significant amound of overhead, optimization is important, but in c++ it is generally 
much faster than python and also pure gdscript implementation. Second, it make me much easiler to 
test the numerical equvialency bettwen these two implementation. In addition, there is no edge case 
check at this stage, I think the gate keeping should be done in a higher level (QuestHandler class in gdscript), 
therefore this object should not be directly used in Godot in production. 
*/

/* 
=> Linear Interpolation Function
Which should be numerically equvilent to np.interp(x, xp, fp)
Parameters
x: array_like (I guess we only need to interpolate one value so double should do the job)
    The x-coordinates at which to evaluate the interpolated values.
xp: 1-D sequence of floats
    The x-coordinates of the data points, must be increasing if argument period is not specified. Otherwise, xp is internally sorted after normalizing the periodic boundaries with xp = xp % period.
fp: 1-D sequence of floats 
    The y-coordinates of the data points, same length as xp.

we use std::lerp which comes from the standard library <cmath.h>
*/

std::pair<int, int> binary_search(std::vector<double>& a, double e, int L, int R) {
    int M = L + (R - L) / 2;

    if (a[M] >= e && a[M-1] <= e) return {M-1, M}; 
    if (a[M] <= e && a[M+1] >= e) return {M, M+1}; 

    if (a[M-1] > e) return binary_search(a, e, L, M-2); 
    if (a[M+1] < e) return binary_search(a, e, M+1, R);

    return {-1, -1};
}

double interp(double x, std::vector<double>& xp, std::vector<double>& fp){

    if (xp.size() != fp.size() || xp.empty()) {
        UtilityFunctions::push_warning("xp and fp must be the same size and non-empty.");
    }

    // some edge cases 
    if (x <= xp.front()) return fp.front();
    if (x >= xp.back()) return fp.back();

    // binary search for intervals NOTE: Now x is bounded b, it's a recursion function, I know it's dangerious, but you can only trust me right now 
    std::pair<int, int> interval = binary_search(xp, x, 0, xp.size()); 

    // Shouldn't failed since we have boundded in front of it but just in case 
    if (interval.first == -1 || interval.second == -1) {
        UtilityFunctions::push_warning("Interpolation failed, Couldn't find correct interval"); 
    }

    // now we have the interval, since 
    double t = (x - xp[interval.first]) / (xp[interval.second] - xp[interval.first]);
    double a = fp[interval.first]; 
    double b = fp[interval.second];

    // now we can call the std::lerp 
    return std::lerp(a, b, t);
}

/*
=> Linear Interpolation Function
*/ 


void Quest::_bind_methods() {
    ClassDB::bind_method(D_METHOD("initialize", "t_guess", "t_guess_sd", "p_threshold", "beta", "delta", "gamma", "grain", "range"), &Quest::initialize, DEFVAL(0.01f), DEFVAL(5));
    ClassDB::bind_method(D_METHOD("quantile", "quantile_order"), &Quest::quantile, DEFVAL(Variant()));
    ClassDB::bind_method(D_METHOD("update", "intensity", "response"), &Quest::update);
    ClassDB::bind_method(D_METHOD("get_pdf"), &Quest::get_pdf);
    ClassDB::bind_method(D_METHOD("get_intensity_history"), &Quest::get_intensity_history);
    ClassDB::bind_method(D_METHOD("get_response_history"), &Quest::get_response_history);
}
Quest::Quest(){}

Quest::~Quest(){}

void Quest::initialize(double t_guess, double t_guess_sd, double p_threshold, double beta, double delta, double gamma, double grain, int range) {
    this->t_guess = t_guess;
    this->t_guess_sd = t_guess_sd;
    this->p_threshold = p_threshold;
    this->beta = beta;
    this->delta = delta;
    this->gamma = gamma;
    this->grain = grain;

    this->dim = (range <= 0) ? 500 : (2 * std::ceil(range / grain / 2.0));
    this->update_pdf = true;
    this->warn_pdf = true;
    this->normalize_pdf = false;

    recompute();
}

void Quest::recompute() {
    // reference: http://courses.washington.edu/matlab1/Lesson_5.html and mostly original source code 
    // this function will be tring to mimic the original function as much as possible 
    
    if (!update_pdf) return; // don't update if user forbid update pdf 
    if (gamma > p_threshold) {
        UtilityFunctions::push_warning(vformat("reducing gamma from %.2f to 0.5", gamma));
        gamma = 0.5;
    }

    i.clear(); for (double val = -dim /2.0; val <= dim/ 2.0; val+= 1.0) i.push_back(val); 
    x.clear(); for (double val: i) x.push_back(val * grain);
    pdf.clear(); for (double val: x) pdf.push_back(std::exp(-0.5 * std::pow(val / t_guess_sd, 2)));

    // normalize the pdf (didn't check if array sum is zero, don't think it's necessary)
    double pdf_sum = std::accumulate(pdf.begin(), pdf.end(), 0.0);
    std::transform(pdf.begin(), pdf.end(), pdf.begin(), [pdf_sum](double val) {return val / pdf_sum;}); // a functional trick

    std::vector<double> i2; for (double val = -dim; val <= dim; val+= 1.0) i2.push_back(val); 
    x2.clear(); for (double val: i2) x2.push_back(val * grain);
    p2.clear(); for (double val: x2) p2.push_back(delta * gamma + (1 - delta) * (1 - (1 - gamma) * std::exp(-std::pow(10, beta * val))));

    std::vector<double> diff(p2.size() -1); // differece between  
    std::vector<int> index; 

    for(size_t ii = 1; ii < p2.size(); ++ii ){
        diff[ii - 1] = p2[ii] - p2[ii-1];
    }

    for(size_t ii = 0; ii < diff.size(); ++ii){
        if (diff[ii] != 0) {
            index.push_back(ii);
        }
    }

    // prepare the xp, fp and x
    std::vector<double> p2_indexed; 
    std::vector<double> x2_indexed; 
    for (int ii: index) {
        p2_indexed.push_back(p2[ii]); 
        x2_indexed.push_back(x2[ii]);
    }

    x_threshold = interp(p_threshold, p2_indexed, x2_indexed); 
    p2.resize(x2.size());
    for (size_t ii = 0; ii < x2.size(); ++ii) {
        double exponent = -std::pow(10, beta * (x2[ii] + x_threshold));
        double exp_value = std::exp(exponent);
        p2[ii] = delta * gamma + (1 - delta) * (1 - (1 - gamma) * exp_value);
    }

    size_t size = p2.size();
    std::vector<double> reversed_p2 = p2;
    std::vector<double> reversed_one_minus_p2(size);
    std::reverse(reversed_p2.begin(), reversed_p2.end());

    for (size_t i = 0; i < size; ++i) {
        reversed_one_minus_p2[i] = 1.0 - p2[size - 1 - i];
    }

    s2 = {reversed_one_minus_p2, reversed_p2};

    double eps = 1e-14; 

    double pL = p2[0]; 
    double pH = p2[p2.size() - 1];
    double pE = pH * std::log(pH + eps) - pL * std::log(pL + eps) + (1-pH+eps)*std::log(1-pH + eps) - (1-pL + eps) * std::log(1-pL+eps);
    pE = 1/(1+std::exp(pE/(pL-pH))); 
    quantile_order = (pE - pL) / (pH - pL); 

    for (size_t ii = 0; ii < intensity_history.size(); ++ii) {
        double inten = std::clamp(intensity_history[ii], static_cast<double>(-1e10), static_cast<double>(1e10));
        std::vector <int> jj; 
        for (double val: i) jj.push_back(static_cast<int>(pdf.size()) + val - std::round((inten - t_guess) / grain) - 1);

        int s2_col = static_cast<int>(s2[0].size()); 
        if (jj[0] < 0 || jj[jj.size() - 1] > s2_col ) {
            if (jj[0] < 0){
                int jj_0 = jj[0];
                std::transform(jj.begin(), jj.end(), jj.begin(), [jj_0](int x) {return x - jj_0;});
            } else {
                int jj_last = jj[jj.size()-1]; 
                std::transform(jj.begin(), jj.end(), jj.begin(),[s2_col, jj_last](int x) {return x + s2_col - jj_last;});
            }

            int response = response_history[ii]; 
            std::vector<double> new_pdf; 
            new_pdf.resize(pdf.size()); 
            for (size_t kk = 0; kk < pdf.size(); ++kk) {
                new_pdf[kk] = pdf[kk] * s2[response][jj[kk]]; // I'm losing my mind, now you know it;s a fucking bad idea to name your variable ii, jj, kk right? Chengran!
            }
            pdf = std::move(new_pdf); 
            /*
                very ambiguous case here 
                the original code: if self.normalizePdf and ii % 100 == 0: 
                since ii is jj in our case, it's an array 
                i'm just gonna assume any now
            */
            if (normalize_pdf && std::any_of(jj.begin(), jj.end(), [](int x) { return x % 100 == 0;})) {
                double pdf_sum = std::accumulate(pdf.begin(), pdf.end(), 0.0);
                std::transform(pdf.begin(), pdf.end(), pdf.begin(), [pdf_sum](double val) {return val / pdf_sum;});
            }

        }

    }

    if (normalize_pdf){
        double pdf_sum = std::accumulate(pdf.begin(), pdf.end(), 0.0);
        std::transform(pdf.begin(), pdf.end(), pdf.begin(), [pdf_sum](double val) {return val / pdf_sum;});
    }

}


double Quest::quantile(const Variant &quantile_order_var) {
    // Step 1: Convert Variant to Optional<double>
    std::optional<double> opt_quantile_order = (quantile_order_var.get_type() != Variant::NIL) 
        ? std::optional<double>(quantile_order_var.operator double()) 
        : std::nullopt;

    // Step 2: Extract quantile order safely
    double quantileOrder = opt_quantile_order.value_or(this->quantile_order);
    if (quantileOrder < 0) quantileOrder = 0.5;

    // not quantileOrder is the one we actually use internally 


    std::vector<double> p(pdf.size());
    std::partial_sum(pdf.begin(), pdf.end(), p.begin());
    std::vector<double> m1p; 
    m1p.push_back(-1);
    m1p.insert(m1p.end(), p.begin(), p.end()); 

    std::vector<double> diff(m1p.size() - 1);
    for (size_t i = 1; i < m1p.size(); ++i) {
        diff[i - 1] = m1p[i] - m1p[i - 1];
    }

    std::vector<int> index;
    for (size_t i = 0; i < diff.size(); ++i) {
        if (diff[i] != 0) {
            index.push_back(i);
        }
    }


    std::vector<double> p_selected;
    std::vector<double> x_selected; 

    for (int idx : index) {
        p_selected.push_back(p[idx]);
        x_selected.push_back(x[idx]);
    }

    double ires = interp(quantileOrder * p[p.size() -1], p_selected, x_selected); 
    
    auto print_vector = [](const std::vector<double> &vec, const String &name) {
        String output = name + String(": [ ");
        for (const auto &val : vec) {
            output += String::num_real(val, 10) + String("\n");  
        }
        output += String("]");
        UtilityFunctions::print(output);
    };

    return t_guess + ires; 
}

void Quest::update(double intensity, int response) {
    if (update_pdf) {
        double inten = std::clamp(intensity, static_cast<double>(-1e10), static_cast<double>(1e10));
        std::vector <int> jj; 
        for (double val: i) jj.push_back(static_cast<int>(pdf.size()) + val - std::round((inten - t_guess) / grain) - 1);

        int s2_col = static_cast<int>(s2[0].size()); 
        if (jj[0] < 0 || jj[jj.size() - 1] > s2_col ) {
            if (jj[0] < 0){
                int jj_0 = jj[0];
                std::transform(jj.begin(), jj.end(), jj.begin(), [jj_0](int x) {return x - jj_0;});
            } else {
                int jj_last = jj[jj.size()-1]; 
                std::transform(jj.begin(), jj.end(), jj.begin(),[s2_col, jj_last](int x) {return x + s2_col - jj_last;});
            }

        }
        
        std::vector<double> new_pdf; 
        new_pdf.resize(pdf.size()); 
        for (size_t kk = 0; kk < pdf.size(); ++kk) {
            new_pdf[kk] = pdf[kk] * s2[response][jj[kk]]; // I'm losing my mind, now you know it's a fucking bad idea to name your variable ii, jj, kk right? Chengran!
        }
        pdf.swap(new_pdf);
        /*
            very ambigous case here 
            the original thing if self.normalizePdf and ii % 100 == 0: 
            since ii is jj in our case, it's an array 
            i'm just gonna assume any now
        */
        if (normalize_pdf) {
            double pdf_sum = std::accumulate(pdf.begin(), pdf.end(), 0.0);
            std::transform(pdf.begin(), pdf.end(), pdf.begin(), [pdf_sum](double val) {return val / pdf_sum;});
        }


    }

    intensity_history.push_back(intensity); 
    response_history.push_back(response); 
}


Array Quest::get_pdf(){
    Array godot_array; 
    for (double val: pdf) godot_array.append(val);
    return godot_array; 
}

Array Quest::get_intensity_history(){
    Array godot_array; 
    for (double val: intensity_history) godot_array.append(val);
    return godot_array;
}


Array Quest::get_response_history(){
    Array godot_array; 
    for (int val: response_history) godot_array.append(val);
    return godot_array;
}

