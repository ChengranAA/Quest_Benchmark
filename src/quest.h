#ifndef QUEST_H
#define QUEST_H

#include <godot_cpp/classes/object.hpp>
#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/variant/utility_functions.hpp> 
#include <godot_cpp/variant/variant.hpp>
#include "godot_cpp/variant/array.hpp" 
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <utility>
#include <optional>

namespace godot {

class Quest: public Object {
    GDCLASS(Quest, Object)

private:

    double t_guess;
    double t_guess_sd;
    double p_threshold;
    double beta;
    double delta;
    double gamma;
    double grain;
    int dim;
    bool update_pdf;
    bool warn_pdf;
    bool normalize_pdf;
    double x_threshold;
    double quantile_order;

    std::vector<double> x;
    std::vector<double> i;
    std::vector<double> pdf;
    std::vector<double> x2;
    std::vector<double> p2;
    std::vector<std::vector<double>> s2;

    std::vector<double> intensity_history;
    std::vector<int> response_history;

    void recompute();

protected:
	static void _bind_methods();

public:
    Quest();
    ~Quest(); 

    void initialize(double t_guess, double t_guess_sd, double p_threshold, double beta, double delta, double gamma, double grain = 0.01f, int range = 5);
    double quantile(const Variant &quantile_order_var);
    void update(double intensity, int response);
    Array get_pdf(); 
    Array get_intensity_history(); 
    Array get_response_history(); 
};

}
#endif