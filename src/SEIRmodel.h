#ifndef SEIR_MODEL_H_
#define SEIR_MODEL_H_

/* NOTE: user-level documentation is provided by corresponding R6 class */

// virtual methods that will be needed for between-farm spread in C++:
class BaseModel
{
public:
  virtual bool is_cpp() const = 0;
  virtual int get_id() const = 0;
  virtual void save() = 0;
  virtual void reset() = 0;
  virtual void set_trans_external(const double trans_external) = 0;
  virtual void update(double d_time) = 0;
  virtual Rcpp::DataFrame get_state() const = 0;
  virtual ~BaseModel() = default;
};

enum class update_type
{
  deterministic,
  stochastic
};

enum class trans_type
{
  frequency,
  density
};

// TODO: clean this up, and specialise for C is array<T, 1> (i.e. redistribute is a no-op)
template<typename T, typename C>
class subcomp;

template<typename C>
class subcomp<double, C>
{
private:
  using T = double;
  void distribute(T total)
  {
    for(int i=0; i<values.size(); ++i)
    {
      values[i] = total / static_cast<double>(values.size());
    }
  }

public:
  C values;

  void set_total(T total)
  {
    distribute(total);
  }

  T get_total() const
  {
    T total = static_cast<T>(0);
    for (const T val : values) total+=val;
    return total;
  }

  void reset()
  {
    distribute(get_total());
  }

};

template<typename C>
class subcomp<int, C>
{
private:
  using T = int;
  void distribute(T total)
  {
    Rcpp::NumericVector probs(values.size(), 1.0/values.size());
    Rcpp::IntegerVector output(values.size());
    R::rmultinom(total, probs.begin(), probs.size(), output.begin());
    for(int i=0; i<values.size(); ++i)
    {
      values[i] = output[i];
    }
  }

public:
  C values;

  void set_total(T total)
  {
    distribute(total);
  }

  T get_total() const
  {
    T total = static_cast<T>(0);
    for (const T val : values) total+=val;
    return total;
  }

  void reset()
  {
    distribute(get_total());
  }

};


template<typename T, typename cE, typename cI, typename cR>
struct state
{
  T S = static_cast<T>(99);
  subcomp<T,cE> Es; // Initialised externally
  subcomp<T,cI> Is; // Initialised externally
  subcomp<T,cR> Rs; // Initialised externally
  T N = static_cast<T>(0);
  double timepoint = 0.0;
};

struct parameters
{
  double beta = 0.05;
  double omega = 0.05;
  double gamma = 0.025;
  double delta = 0.005;
  double vacc = 0.0;//0.001;
  double repl = 0.0;//0.0001;
  double cull = 0.0;//0.002;
  double trans_external = 0.0;
  trans_type transmission_type = trans_type::frequency;
};


template<bool t_debug>
void assertr(bool cond, std::string_view msg);
template <>
void assertr<false>(bool cond, std::string_view msg) { }

// TODO: decompose the specialisation so there are fewer combinations???
// templated class for specialisation (stochastic/deterministic, sir/seir/numE=0vs3vsX??, debug)
template<update_type t_update_type, int t_numE, int t_numI, int t_numR, bool t_debug>
class SEIRmodel : public BaseModel
{
private:
  // TODO: some way to make this constexpr conditional on t_fixedE???
  using t_numE_type = const int;
  // typedef std::conditional<t_fixedE, constexpr int, const int>::type t_numE_type;
  t_numE_type m_numE;
  t_numE_type m_numI;
  t_numE_type m_numR;

  static constexpr bool t_fixedE = t_numE>=0;
  static constexpr bool t_hasE = t_numE!=0;
  static constexpr bool t_fixedI = t_numI>=0;
  static constexpr bool t_hasI = t_numI!=0;
  static constexpr bool t_fixedR = t_numR>=0;
  static constexpr bool t_hasR = t_numR!=0;

  std::string m_group_name;
  bool m_has_name;

  trans_type m_trans_type = trans_type::frequency;

  parameters m_pars;
  parameters m_saved_pars;

  typedef typename std::conditional<
      t_update_type == update_type::deterministic, double, int
    >::type t_state_type;
  typedef typename std::conditional<
      t_fixedE, std::array<t_state_type, (t_numE<0 ) ? 0 : t_numE>, std::vector<t_state_type>
    >::type t_Etype;
  typedef typename std::conditional<
      t_fixedI, std::array<t_state_type, (t_numI<0 ) ? 0 : t_numI>, std::vector<t_state_type>
    >::type t_Itype;
  typedef typename std::conditional<
      t_fixedR, std::array<t_state_type, (t_numR<0 ) ? 0 : t_numR>, std::vector<t_state_type>
    >::type t_Rtype;
  state<t_state_type, t_Etype, t_Itype, t_Rtype> m_state;
  state<t_state_type, t_Etype, t_Itype, t_Rtype> m_saved_state;

  int m_iteration = 0;

public:
  SEIRmodel(const int numE, const int numI, const int numR, const Rcpp::StringVector group_name) :
    m_numE(t_fixedE ? t_numE : numE),
    m_numI(t_fixedI ? t_numI : numI),
    m_numR(t_fixedR ? t_numR : numR),
    m_group_name(Rcpp::as< std::string > (group_name(0L))),
    m_has_name(!Rcpp::is_na(group_name)[0L])
  {

    if constexpr (t_debug){
      if (t_fixedE && (t_numE != numE)) Rcpp::stop("Invalid numE");
      if (t_fixedI && (t_numI != numI)) Rcpp::stop("Invalid numI");
      if (t_fixedR && (t_numR != numR)) Rcpp::stop("Invalid numR");
      if (group_name.size()!=1L) Rcpp::stop("Invalid group_name");
    }

    if constexpr (!t_hasE) {
      m_pars.omega = 0.0;
    } else {
      if constexpr (!t_fixedE){
        m_state.Es.values.resize(m_numE);
      }
      m_state.Es.set_total(static_cast<t_state_type>(0));
    }

    // We always have at least 1 I!
    if constexpr (!t_fixedI){
      m_state.Is.values.resize(m_numI);
    }
    m_state.Is.set_total(static_cast<t_state_type>(1));

    if constexpr (!t_hasR) {
      m_pars.delta = 0.0;
      m_pars.vacc = 0.0;
    } else {
      if constexpr (!t_fixedR){
        m_state.Rs.values.resize(m_numR);
      }
      m_state.Rs.set_total(static_cast<t_state_type>(0));
    }

    reset_N();

    // Run save method:
    save();

    // Test reset method:
    reset();
  }

  // Overridden virtual methods:
  bool is_cpp() const noexcept override
  {
    return true;
  }

  int get_id() const override
  {
    // TODO: generate ID in constructor (maybe in base class?)
    return 0;
  }

  void save() override
  {
    m_saved_state = m_state;
    m_saved_pars = m_pars;
  }

  void reset() override
  {
    m_state = m_saved_state;
    // For compatibility with R6:
    if constexpr (t_hasE) m_state.Es.reset();
    m_state.Is.reset();
    if constexpr (t_hasR) m_state.Rs.reset();
    m_pars = m_saved_pars;
  }

  /*
  void set_trans_external(double trans_external) override
  {

  }*/

  double transmission_rate()
  {
    double val = 0.0;
    if(m_pars.transmission_type == trans_type::frequency){
      val = m_pars.beta * m_state.Is.get_total() / m_state.N;
    }else if(m_pars.transmission_type == trans_type::density){
      val = m_pars.beta * m_state.Is.get_total();
    }else{
      Rcpp::stop("Unrecognised m_pars.transmission_type");
    }
    return val + get_trans_external();
  }

  /*
  void update()
  {
    update(m_d_time);
  }*/

  void update(double d_time) override
  {
    if(d_time <= 0.0) Rcpp::stop("Invalid d_time");

    // Keep the old version here:
    state old_state = m_state;

    std::array<double, 2> srates = { transmission_rate(), m_pars.vacc };
    std::array<t_state_type, 2> leave_S = apply_rates<2>(old_state.S, srates, d_time);
    m_state.S -= (leave_S[0] + leave_S[1]);
    if constexpr (t_hasE) {
      m_state.Es.values[0] += leave_S[0];
    }else{
      m_state.Is.values[0] += leave_S[0];
    }
    if constexpr (t_hasR) {
      m_state.Rs.values[0] += leave_S[1];
    }

    if constexpr (t_hasE) {
      t_state_type carry = static_cast<t_state_type>(0);
      t_state_type toS = static_cast<t_state_type>(0);
      std::array<double, 2> rates = { m_pars.omega*m_numE, m_pars.repl };
      for(int i=0; i<(t_numE>0 ? t_numE : m_numE); ++i){
        std::array<t_state_type, 2> leave = apply_rates<2>(old_state.Es.values[i], rates, d_time);
        m_state.Es.values[i] += carry - leave[0] - leave[1];
        carry = leave[0];
        toS += leave[1];
      }
      m_state.S += toS;
      m_state.Is.values[0] += carry;
    }

    // We always have I!
    {
      t_state_type carry = static_cast<t_state_type>(0);
      t_state_type toS = static_cast<t_state_type>(0);
      std::array<double, 2> rates = { m_pars.gamma*m_numI, m_pars.repl + m_pars.cull };
      for(int i=0; i<(t_numI>0 ? t_numI : m_numI); ++i){
        std::array<t_state_type, 2> leave = apply_rates<2>(old_state.Is.values[i], rates, d_time);
        m_state.Is.values[i] += carry - leave[0] - leave[1];
        carry = leave[0];
        toS += leave[1];
      }
      m_state.S += toS;
      if constexpr (t_hasR) {
        m_state.Rs.values[0] += carry;
      } else {
        m_state.S += carry;
      }
    }

    if constexpr (t_hasR) {
      t_state_type carry = static_cast<t_state_type>(0);
      t_state_type toS = static_cast<t_state_type>(0);
      std::array<double, 2> rates = { m_pars.delta*m_numR, m_pars.repl };
      for(int i=0; i<(t_numR>0 ? t_numR : m_numR); ++i){
        std::array<t_state_type, 2> leave = apply_rates<2>(old_state.Rs.values[i], rates, d_time);
        m_state.Rs.values[i] += carry - leave[0] - leave[1];
        carry = leave[0];
        toS += leave[1];
      }
      m_state.S += toS;
      m_state.S += carry;
    }

    m_state.timepoint += d_time;
    check_state();
  }

  Rcpp::DataFrame get_state() const override
  {
    Rcpp::DataFrame rv = Rcpp::DataFrame::create(
      Rcpp::Named("Time") = m_state.timepoint,
      Rcpp::Named("S") = m_state.S,
      Rcpp::Named("E") = m_state.Es.get_total(),
      Rcpp::Named("I") = m_state.Is.get_total(),
      Rcpp::Named("R") = m_state.Rs.get_total()
    );

    if(m_has_name)
    {
      Rcpp::String group_name = m_group_name;
      rv.push_front(group_name, "Group");
    }
    if constexpr (!t_hasR) {
      rv.erase(4);
    }
    if constexpr (!t_hasE) {
      rv.erase(2);
    }

    return Rcpp::as<Rcpp::DataFrame>(rv);
  }

  // New methods:
  template<int t_size>
  std::array<t_state_type, t_size> apply_rates(const t_state_type compartment, const std::array<double, t_size>& rates, const double d_time)
  {
    double sumrates = 0.0;
    for(double r : rates) sumrates+=r;
    const double leave = 1.0 - std::exp(-sumrates*d_time);

    Rcpp::NumericVector probs(t_size+1);
    double sumprob = 0.0;
    for(int i=0; i<t_size; ++i){
      probs[i+1] = leave==0.0 ? 0.0 : (leave * rates[i] / sumrates);
      sumprob += probs[i+1];
    }
    probs[0] = 1.0 - sumprob;

    std::array<t_state_type, t_size> rv;
    if constexpr (t_update_type==update_type::deterministic)
    {
      for(int i=0; i<t_size; ++i){
        rv[i] = compartment * probs[i+1];
      }
    }
    else if constexpr (t_update_type==update_type::stochastic)
    {
      Rcpp::IntegerVector output(t_size+1);
      R::rmultinom(compartment, probs.begin(), t_size+1, output.begin());
      for(int i=0; i<t_size; ++i){
        rv[i] = output[i+1];
      }
    }
    else{
      Rcpp::stop("Unhandled t_update_type");
    }
    return rv;
  }


  void reset_N()
  {
    m_state.N = m_state.S+m_state.Es.get_total()+m_state.Is.get_total()+m_state.Rs.get_total();
    check_state();
  }

  void set_S(t_state_type val)
  {
    m_state.S = val;
    reset_N();
  }
  t_state_type get_S() const
  {
    return m_state.S;
  }

  void set_E(t_state_type val)
  {
    if constexpr (!t_hasE) Rcpp::stop("Model contains no E compartments!");
    m_state.Es.set_total(val);
    reset_N();
  }
  t_state_type get_E() const
  {
    return m_state.Es.get_total();
  }

  void set_I(t_state_type val)
  {
    m_state.Is.set_total(val);
    reset_N();
  }
  t_state_type get_I() const
  {
    return m_state.Is.get_total();
  }

  void set_R(t_state_type val)
  {
    if constexpr (!t_hasR) Rcpp::stop("Model contains no R compartments!");
    m_state.Rs.set_total(val);
    reset_N();
  }
  t_state_type get_R() const
  {
    return m_state.Rs.get_total();;
  }

  t_state_type get_N() const
  {
    return m_state.N;
  }

  void set_beta(double val)
  {
    m_pars.beta = val;
  }
  double get_beta() const
  {
    return m_pars.beta;
  }

  void set_omega(double val)
  {
    if constexpr(!t_hasE) Rcpp::stop("Unable to change omega:  no E state!");
    m_pars.omega = val;
  }
  double get_omega() const
  {
    return m_pars.omega;
  }

  void set_gamma(double val)
  {
    m_pars.gamma = val;
  }
  double get_gamma() const
  {
    return m_pars.gamma;
  }

  void set_delta(double val)
  {
    if constexpr(!t_hasR) Rcpp::stop("Unable to change delta:  no R state!");
    m_pars.delta = val;
  }
  double get_delta() const
  {
    return m_pars.delta;
  }

  void set_repl(double val)
  {
    m_pars.repl = val;
  }
  double get_repl() const
  {
    return m_pars.repl;
  }

  void set_cull(double val)
  {
    m_pars.cull = val;
  }
  double get_cull() const
  {
    return m_pars.cull;
  }

  void set_vacc(double val)
  {
    if constexpr(!t_hasR) Rcpp::stop("Unable to change vacc:  no R state!");
    m_pars.vacc = val;
  }
  double get_vacc() const
  {
    return m_pars.beta;
  }

  double get_time() const
  {
    return m_state.timepoint;
  }

  void set_trans_external(double trans_external) override
  {
    m_pars.trans_external = trans_external;
  }
  double get_trans_external() const
  {
    return m_pars.trans_external;
  }

  void set_trans_type(Rcpp::String val)
  {
    if(val=="frequency"){
      m_pars.transmission_type = trans_type::frequency;
    }else if(val=="density"){
      m_pars.transmission_type = trans_type::density;
    }else{
      Rcpp::stop("Unrecognised transmission_type");
    }
  }
  Rcpp::String get_trans_type() const
  {
    Rcpp::String val;
    if(m_pars.transmission_type == trans_type::frequency){
      val = "frequency";
    }else if(m_pars.transmission_type == trans_type::density){
      val = "density";
    }else{
      Rcpp::stop("Unrecognised m_pars.transmission_type");
    }
    return val;
  }

  void check_state() const
  {
    if constexpr (t_debug)
    {
      t_state_type newN = m_state.S+m_state.Es.get_total()+m_state.Is.get_total()+m_state.Rs.get_total();
      if constexpr (t_update_type==update_type::deterministic){
        if(std::abs(newN - m_state.N) > 1e-6) Rcpp::stop("Error in update: N has changed!");
      }else{
        if(newN != m_state.N) Rcpp::stop("Error in update: N has changed!");
      }
      // TODO: check parameter values and state>0 etc
    }
  }

  void show() const
  {

    Rcpp::Rcout << "An SEIR model with ";
    if (m_has_name) Rcpp::Rcout << "identifier/name '" << m_group_name << "' and ";
    Rcpp::Rcout << "the following properties:\n\t";
    Rcpp::Rcout << "S/E/I/R (N) = " << m_state.S << "/" << m_state.Es.get_total() << "/" << m_state.Is.get_total() << "/" << m_state.Rs.get_total() << " (" << m_state.N << ")\n\t";
    Rcpp::Rcout << "beta/omega/gamma/delta = " << m_pars.beta << "/" << m_pars.omega << "/" << m_pars.gamma << "/" << m_pars.delta << "\n\t";
    Rcpp::Rcout << "vacc/repl/cull = " << m_pars.vacc << "/" << m_pars.repl << "/" << m_pars.cull << "\n\t";
    Rcpp::Rcout << "E compartments = " << m_numE << "\n\t";
    Rcpp::Rcout << "I compartments = " << m_numI << "\n\t";
    Rcpp::Rcout << "R compartments = " << m_numR << "\n\t";
    Rcpp::Rcout << "external transmission = " << m_pars.trans_external << "\n\t";
    if constexpr (t_update_type==update_type::deterministic){
      Rcpp::Rcout << "update type = " << "deterministic" << "\n\t";
    }else{
      Rcpp::Rcout << "update type = " << "stochastic" << "\n\t";
    }
    if(m_pars.transmission_type==trans_type::frequency){
      Rcpp::Rcout << "transmission type = " << "frequency" << "\n\t";
    }else{
      Rcpp::Rcout << "transmission type = " << "density" << "\n\t";
    }
    Rcpp::Rcout << "current time = " << m_state.timepoint << "\n";

    /*
    if constexpr (t_debug){
      Rcpp::Rcout << "\t[C++ implementation; DEBUG mode]\n";
    }else{
      Rcpp::Rcout << "\t[C++ implementation]\n";
    }*/

  }

  Rcpp::DataFrame run(double add_time, double d_time)
  {

    const bool include_current = m_state.timepoint==0.0;
    const int addi = static_cast<int>(include_current);
    int rows = static_cast<int>(add_time/d_time) + addi;

    Rcpp::NumericVector Time(rows);
    typedef typename std::conditional<
      t_update_type == update_type::deterministic, Rcpp::NumericVector, Rcpp::IntegerVector
    >::type VT;
    VT S(rows);
    VT E(rows);
    VT I(rows);
    VT R(rows);

    if(include_current){
      constexpr int i=0;
      Time[i] = get_time();
      S[i] = get_S();
      E[i] = get_E();
      I[i] = get_I();
      R[i] = get_R();
    }

    for(int i=addi; i<rows; ++i)
    {
      update(d_time);
      Time[i] = get_time();
      S[i] = get_S();
      E[i] = get_E();
      I[i] = get_I();
      R[i] = get_R();
    }

    Rcpp::DataFrame rv = Rcpp::DataFrame::create(
      Rcpp::Named("Time") = Time,
      Rcpp::Named("S") = S,
      Rcpp::Named("E") = E,
      Rcpp::Named("I") = I,
      Rcpp::Named("R") = R
    );

    if(m_has_name)
    {
      Rcpp::String group_name = m_group_name;
      rv.push_front(group_name, "Group");
    }
    if constexpr (!t_hasR) {
      rv.erase(4);
    }
    if constexpr (!t_hasE) {
      rv.erase(2);
    }

    return Rcpp::as<Rcpp::DataFrame>(rv);
  }

};


#endif // SEIR_MODEL_H_
