#include <iomanip>
#include <sstream>
#include <cmath>

#include "parameter_study.h"

Range::Range():
    start(0.0),
    end(0.0),
    num_steps(0)
{ 
    LOG_ERROR("Range is not initialized witha ny meaningful value.");
}


Range::Range(double start, double end, size_t num_steps):
    start(start),
    end(end),
    num_steps(num_steps)
{
    if (start > end)
    {
        LOG_ERROR("Start value (" + std::to_string(start) + ") is greater than end value (" + std::to_string(end) + ")");
        this->end = this->start;
    }
}


std::string Range::to_array() const
{
    std::stringstream ss;
    ss << "{";

    size_t i_max;
    if (this->num_steps < 10)
        i_max = this->num_steps;
    else
        i_max = 7;
    for (size_t i = 0; i < i_max; ++i)
    {
        ss << this->operator[](i);
        if (i < i_max - 1)
        {
            ss << ", ";
        }
    }
    if (this->num_steps >= 10)
    {
        ss << ", ..., ";
        ss << this->operator[](this->num_steps - 2);
        ss << ", ";
        ss << this->operator[](this->num_steps - 1);
    }
    ss << "}";
    return ss.str();
}


std::ostream &operator<<(std::ostream &os, const Range &range)
{
    os << range.to_string();
    return os;
}


size_t Range::get_num_steps() const
{
    return this->num_steps;
}


Const::Const(double value):
    Range(value, value, 1)
{}


double Const::operator[](size_t i) const
{
    (void)i;
    return this->start;
}


std::string Const::to_string() const
{
    return "Const(" + std::to_string(this->start) + ")";
}


LinearRange::LinearRange(double start, double end, size_t num_steps):
    Range(start, end, std::max(num_steps, (size_t)1))
{}


double LinearRange::operator[](size_t i) const
{
    if (this->num_steps <= 1)
    {
        return this->start;
    }
    if (i >= this->num_steps)
    {
        return this->end;
    }

    return this->start + i * (this->end - this->start) / (this->num_steps - 1);
}


std::string LinearRange::to_string() const
{
    return "LinearRange(" + std::to_string(this->start) + ", " + std::to_string(this->end) + ", " + std::to_string(this->num_steps) + ")";
}


PowRange::PowRange(double start, double end, size_t num_steps, double base):
    Range(start, end, std::max(num_steps, (size_t)1))
{
    if (base <= 0.0)
    {
        LOG_ERROR("The value of base (" + std::to_string(base) + ") must be positive. Otherwise the formula f(x) = a * x^base + c breaks.");
        this->start = 1.0;
    }
    this->base = base;
    this->a = (end - start) / std::pow(num_steps - 1, base);
    this->c = start;
}


double PowRange::operator[](size_t i) const
{
    if (this->num_steps <= 1)
    {
        return this->start;
    }
    if (i >= this->num_steps)
    {
        return this->end;
    }

    return this->a * std::pow(i, this->base) + this->c;
}


std::string PowRange::to_string() const
{
    return "PowRange(" + std::to_string(this->start) + ", " + std::to_string(this->end) + ", " + std::to_string(this->num_steps) + ", " + std::to_string(this->base) + ")";
}


std::unique_ptr<Range> get_unique_ptr(const ParameterStudy::AnyRange range)
{
    if (std::holds_alternative<Const>(range))
    {
        return std::make_unique<Const>(std::get<Const>(range));
    }
    else if (std::holds_alternative<LinearRange>(range))
    {
        return std::make_unique<LinearRange>(std::get<LinearRange>(range));
    }
    else if (std::holds_alternative<PowRange>(range))
    {
        return std::make_unique<PowRange>(std::get<PowRange>(range));
    }
    return nullptr;
}

ParameterStudy::ParameterStudy(const ParameterStudy::Builder &builder):
    combination_ID(0),
    mechanism(builder.mechanism),
    R_E(get_unique_ptr(builder.R_E)),
    species(builder.species),
    fractions(builder.fractions),
    P_amb(get_unique_ptr(builder.P_amb)),
    T_inf(get_unique_ptr(builder.T_inf)),
    alfa_M(get_unique_ptr(builder.alfa_M)),
    P_v(get_unique_ptr(builder.P_v)),
    mu_L(get_unique_ptr(builder.mu_L)),
    rho_L(get_unique_ptr(builder.rho_L)),
    c_L(get_unique_ptr(builder.c_L)),
    surfactant(get_unique_ptr(builder.surfactant)),
    enable_heat_transfer(builder.enable_heat_transfer),
    enable_evaporation(builder.enable_evaporation),
    enable_reactions(builder.enable_reactions),
    enable_dissipated_energy(builder.enable_dissipated_energy),
    target_specie(builder.target_specie),
    excitation_type(builder.excitation_type)
{
    this->excitation_params.reserve(builder.excitation_params.size());
    for (const auto &excitation_param : builder.excitation_params)
    {
        this->excitation_params.push_back(get_unique_ptr(excitation_param));
    }
}


std::string ParameterStudy::to_string(const bool with_code) const
{
    const Parameters* par = Parameters::get_parameters(this->mechanism);
    if (par == nullptr)
    {
        return "";
    }
    std::stringstream ss;
    const size_t strw = 28;
    const size_t rangew = 48;
    auto format_field_name = [](std::ostream& os) -> std::ostream& {
        return os << "    " << std::setw(strw);
    };
    auto format_range = [](Range *range, const bool last = false) -> std::string {
        std::stringstream ss;
        ss << std::left << std::setw(rangew) << (range->to_string() + (last ? " " : ",")) << "    // " << range->to_array();
        return ss.str();
    };
    auto species_to_string = [&with_code](const std::string name) -> std::string {
        if (with_code) return "\"" + name + "\"";
        else return name;
    };
    std::vector<std::string> species_names;
    species_names.reserve(this->species.size());
    for (auto &species_name : this->species)
    {
        species_names.push_back(species_to_string(species_name));
    }

    if (with_code) ss << "ParameterStudy::Builder{\n";
    ss << std::left;

    ss << format_field_name << ".mechanism" << " = " << ((with_code ? "Parameters::mechanism::" : "") + par->model) << ",\n";
    ss << format_field_name << ".R_E" << " = " << format_range(this->R_E.get()) << "\n";
    ss << format_field_name << ".species" << " = " << ::to_string((std::string*)species_names.data(), species_names.size()) << ",\n";
    ss << format_field_name << ".fractions" << " = " << ::to_string((double*)this->fractions.data(), this->fractions.size()) << ",\n";
    ss << format_field_name << ".P_amb" << " = " << format_range(this->P_amb.get()) << "\n";
    ss << format_field_name << ".T_inf" << " = " << format_range(this->T_inf.get()) << "\n";
    ss << format_field_name << ".alfa_M" << " = " << format_range(this->alfa_M.get()) << "\n";
    ss << format_field_name << ".P_v" << " = " << format_range(this->P_v.get()) << "\n";
    ss << format_field_name << ".mu_L" << " = " << format_range(this->mu_L.get()) << "\n";
    ss << format_field_name << ".rho_L" << " = " << format_range(this->rho_L.get()) << "\n";
    ss << format_field_name << ".c_L" << " = " << format_range(this->c_L.get()) << "\n";
    ss << format_field_name << ".surfactant" << " = " << format_range(this->surfactant.get()) << "\n";
    ss << format_field_name << ".enable_heat_transfer" << " = " << std::boolalpha << this->enable_heat_transfer << ",\n";
    ss << format_field_name << ".enable_evaporation" << " = " << std::boolalpha << this->enable_evaporation << ",\n";
    ss << format_field_name << ".enable_reactions" << " = " << std::boolalpha << this->enable_reactions << ",\n";
    ss << format_field_name << ".enable_dissipated_energy" << " = " << std::boolalpha << this->enable_dissipated_energy << ",\n";
    ss << format_field_name << ".target_specie" << " = " << this->target_specie << ",\n";
    ss << format_field_name << ".excitation_params" << " = {";
    for (size_t i = 0; i < this->excitation_params.size(); ++i)
    {
        if (i == 0)
            ss << "\n";
        ss << "        " << excitation_params[i]->to_string();
        if (i < this->excitation_params.size() - 1)
            ss << ",\n";
        else
            ss << "\n";
    }

    ss << "    },\n";
    ss << format_field_name << ".excitation_type"            << " = " << (with_code ? "Parameters::excitation::" : "") << Parameters::excitation_names[this->excitation_type] << "\n";

    if (with_code) ss << "}";
    ss << std::right;
    return ss.str();
}


std::ostream &operator<<(std::ostream &os, const ParameterStudy &ps)
{
    os << ps.to_string(true);
    return os;
}


size_t ParameterStudy::get_total_combination_count() const
{
    size_t count = 1;
    count *= this->R_E->get_num_steps();
    count *= this->P_amb->get_num_steps();
    count *= this->T_inf->get_num_steps();
    count *= this->alfa_M->get_num_steps();
    count *= this->P_v->get_num_steps();
    count *= this->mu_L->get_num_steps();
    count *= this->rho_L->get_num_steps();
    count *= this->c_L->get_num_steps();
    count *= this->surfactant->get_num_steps();
    for (const auto &excitation_param : this->excitation_params)
    {
        count *= excitation_param->get_num_steps();
    }

    return count;
}


size_t ParameterStudy::get_next_combination_ID() const
{
    return this->combination_ID.load();
}


std::pair<is_success, ControlParameters> ParameterStudy::get_next_combination()
{
    const size_t combination_ID = this->combination_ID.fetch_add(1, std::memory_order_relaxed);    
    size_t prod = 1;

    auto get_value = [&combination_ID, &prod](const std::unique_ptr<Range> &range) -> double {
        double value = range->operator[](combination_ID / prod % range->get_num_steps());
        prod *= range->get_num_steps();
        return value;
    };

    ControlParameters cpar{ControlParameters::Builder{
        .ID = combination_ID,
        .mechanism = this->mechanism,
        .R_E = get_value(this->R_E),
        .species = this->species,
        .fractions = this->fractions,
        .P_amb = get_value(this->P_amb),
        .T_inf = get_value(this->T_inf),
        .alfa_M = get_value(this->alfa_M),
        .P_v = get_value(this->P_v),
        .mu_L = get_value(this->mu_L),
        .rho_L = get_value(this->rho_L),
        .c_L = get_value(this->c_L),
        .surfactant = get_value(this->surfactant),
        .enable_heat_transfer = this->enable_heat_transfer,
        .enable_evaporation = this->enable_evaporation,
        .enable_reactions = this->enable_reactions,
        .enable_dissipated_energy = this->enable_dissipated_energy,
        .target_specie = this->target_specie
    }};
    
    for (size_t i = 0; i < this->excitation_params.size(); ++i)
    {
        cpar.excitation_params[i] = get_value(this->excitation_params[i]);
    }

    is_success success = (combination_ID / prod % 2) == 0;
    if (!success)
    {
        std::string message = "combination_ID (" + std::to_string(combination_ID) + ") is out of range for current parameter study with a max combination count of " + std::to_string(this->get_total_combination_count());
        cpar.error_ID = LOG_ERROR(message, combination_ID);
    }

    return {success, cpar};
}