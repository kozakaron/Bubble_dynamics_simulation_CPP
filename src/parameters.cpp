#include "common.h"
#include "parameters.h"

#include "chemkin_ar_he.h"
#include "chemkin_kaust2023_n2.h"
#include "chemkin_otomo2018_without_o.h"
#include "chemkin_otomo2018.h"

const Parameters Parameters::chemkin_ar_he_params                = Parameters(chemkin_ar_he_struct());
const Parameters Parameters::chemkin_kaust2023_n2_params         = Parameters(chemkin_kaust2023_n2_struct());
const Parameters Parameters::chemkin_otomo2018_without_o_params  = Parameters(chemkin_otomo2018_without_o_struct());
const Parameters Parameters::chemkin_otomo2018_params            = Parameters(chemkin_otomo2018_struct());


#define COPY_ARRAY(type, name, size) { \
    type* temp = new type[size]; \
    this->name = (const type*)temp; \
    std::copy((type*)(T::name), (type*)(T::name) + size, temp); \
}


template <typename T>
Parameters::Parameters(T dummy):
    model(T::model),
    input_file(T::input_file),
    num_elements(T::num_elements),
    num_species(T::num_species),
    index_of_water(T::index_of_water),
    invalid_index(T::invalid_index),
    _elements(),
    _species(),
    species_names(),
    NASA_order(T::NASA_order),
    num_reactions(T::num_reactions),
    num_max_specie_per_reaction(T::num_max_specie_per_reaction),
    num_third_bodies(T::num_third_bodies),
    num_irreversible(T::num_irreversible),
    num_pressure_dependent(T::num_pressure_dependent),
    num_lindemann(T::num_lindemann),
    num_troe(T::num_troe),
    num_sri(T::num_sri),
    num_plog(T::num_plog),
    num_plog_levels(T::num_plog_levels)
{
    (void)dummy;
        
    COPY_ARRAY(double, W, T::num_species);
    COPY_ARRAY(double, lambdas, T::num_species);
    COPY_ARRAY(double, temp_range, T::num_species*3);
    COPY_ARRAY(double, a_low, T::num_species*(T::NASA_order+2));
    COPY_ARRAY(double, a_high, T::num_species*(T::NASA_order+2));
    COPY_ARRAY(double, A, T::num_reactions);
    COPY_ARRAY(double, b, T::num_reactions);
    COPY_ARRAY(double, E, T::num_reactions);
    COPY_ARRAY(index_t, nu_indexes, T::num_reactions*T::num_max_specie_per_reaction);
    for(index_t i = 0; i < T::num_elements; i++)
    {
        _elements[T::elements[i].first] = T::elements[i].second;
    }
    for(index_t i = 0; i < T::num_species; i++)
    {
        _species[T::species[i].first] = T::species[i].second;
        species_names.push_back(T::species[i].first);
    }
    COPY_ARRAY(stoich_t, nu, T::num_reactions*T::num_max_specie_per_reaction*3);
    COPY_ARRAY(index_t, third_body_indexes, T::num_third_bodies);
    COPY_ARRAY(bool, is_pressure_dependent, T::num_third_bodies);
    COPY_ARRAY(double, alfa, T::num_third_bodies*T::num_species);
    COPY_ARRAY(index_t, irreversible_indexes, T::num_irreversible);
    COPY_ARRAY(index_t, pressure_dependent_indexes, T::num_pressure_dependent);
    COPY_ARRAY(Parameters::reac_type, pressure_dependent_reac_types, T::num_pressure_dependent);
    COPY_ARRAY(index_t, is_third_body_indexes, T::num_pressure_dependent);
    COPY_ARRAY(double, reac_const, T::num_pressure_dependent*3);
    COPY_ARRAY(double, troe, T::num_troe*4);
    COPY_ARRAY(double, sri, T::num_sri*5);
    COPY_ARRAY(index_t, plog_indexes, T::num_plog);
    if (T::num_plog > 0)
    {
        COPY_ARRAY(index_t, plog_seperators, T::num_plog+1);
    } else {
        COPY_ARRAY(index_t, plog_seperators, 2);
    }
   COPY_ARRAY(double, plog, T::num_plog_levels*4);
}


Parameters::~Parameters()
{
    if (W != nullptr) delete[] W;
    if (lambdas != nullptr) delete[] lambdas;
    if (temp_range != nullptr) delete[] temp_range;
    if (a_low != nullptr) delete[] a_low;
    if (a_high != nullptr) delete[] a_high;
    if (A != nullptr) delete[] A;
    if (b != nullptr) delete[] b;
    if (E != nullptr) delete[] E;
    if (nu_indexes != nullptr) delete[] nu_indexes;
    //if (nu_forward != nullptr) delete[] nu_forward;
    //if (nu_backward != nullptr) delete[] nu_backward;
    if (nu != nullptr) delete[] nu;
    if (third_body_indexes != nullptr) delete[] third_body_indexes;
    if (is_pressure_dependent != nullptr) delete[] is_pressure_dependent;
    if (alfa != nullptr) delete[] alfa;
    if (irreversible_indexes != nullptr) delete[] irreversible_indexes;
    if (pressure_dependent_indexes != nullptr) delete[] pressure_dependent_indexes;
    if (pressure_dependent_reac_types != nullptr) delete[] pressure_dependent_reac_types;
    if (is_third_body_indexes != nullptr) delete[] is_third_body_indexes;
    if (reac_const != nullptr) delete[] reac_const;
    if (troe != nullptr) delete[] troe;
    if (sri != nullptr) delete[] sri;
    if (plog_indexes != nullptr) delete[] plog_indexes;
    if (plog_seperators != nullptr) delete[] plog_seperators;
    if (plog != nullptr) delete[] plog;
}

const Parameters *Parameters::get_parameters(const Parameters::mechanism mech)
{
    switch (mech)
    {
        case mechanism::chemkin_ar_he:
            return &Parameters::chemkin_ar_he_params;
        case mechanism::chemkin_kaust2023_n2:
            return &Parameters::chemkin_kaust2023_n2_params;
        case mechanism::chemkin_otomo2018_without_o:
            return &Parameters::chemkin_otomo2018_without_o_params;
        case mechanism::chemkin_otomo2018:
            return &Parameters::chemkin_otomo2018_params;
        default:
            LOG_ERROR("Unknown mechanisms", 0);
            return nullptr;
    }
}

index_t Parameters::get_element(const string &name) const
{
    auto it = _elements.find(name);
    if (it == _elements.end())
    {
        LOG_ERROR("Element \"" + name + "\" not in " + this->model, 0);
        return invalid_index;
    }
    return it->second;
}

index_t Parameters::get_species(const string &name) const
{
    auto it = _species.find(name);
    if (it == _species.end())
    {
        LOG_ERROR("Species \"" + name + "\" not in " + this->model, 0);
        return invalid_index;
    }
    return it->second;
}