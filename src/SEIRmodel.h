#ifndef SEIR_MODEL_H_
#define SEIR_MODEL_H_

/* NOTE: user-level documentation is provided by corresponding R6 class */

// virtual methods that will be needed for between-farm spread in C++:
class BaseModel
{
public:
  virtual void get_id() const = 0;
  virtual void save() = 0;
  virtual void reset() = 0;
  virtual void set_trans_external(const double trans_external) = 0;
  virtual void update(double d_time) = 0;
  virtual Rcpp::DataFrame get_state() const = 0;
  virtual ~BaseModel() = default;
};

template<bool t_debug>
class SEIRmodel : public BaseModel
{
public:
  SEIRmodel(const int numE, const double d_time, const Rcpp::CharacterVector group_name)
  {
    if constexpr (t_debug){
      Rcpp::Rcout << "DEBUG\n";
    }
    else
    {
      Rcpp::Rcout << "NORMAL\n";
    }
  }

  // Overridden virtual methods:
  void get_id() const override
  {

  }

  void save() override
  {

  }

  void reset() override
  {

  }

  void set_trans_external(double trans_external) override
  {

  }

  void update(double d_time) override
  {

  }

  Rcpp::DataFrame get_state() const override
  {
    return Rcpp::DataFrame();
  }

  // New methods:
  void set_S(int S=1)
  {

  }

  void show() const
  {
    if constexpr (t_debug){
      Rcpp::Rcout << "show DEBUG\n";
    }
    else
    {
      Rcpp::Rcout << "show NORMAL\n";
    }

  }

};


#endif // SEIR_MODEL_H_
