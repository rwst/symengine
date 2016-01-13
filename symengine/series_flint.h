#ifndef SYMENGINE_SERIES_FLINT_H
#define SYMENGINE_SERIES_FLINT_H

#include <symengine/series.h>
#include <symengine/rational.h>
#include <symengine/expression.h>

#ifdef HAVE_SYMENGINE_FLINT
#include <flint/fmpq_polyxx.h>

namespace SymEngine {

using fp_t = flint::fmpq_polyxx;
// Univariate Rational Coefficient Power SeriesBase using Flint
class URatPSeriesFlint : public SeriesBase<fp_t, flint::fmpqxx, URatPSeriesFlint, false> {
public:
    URatPSeriesFlint(const fp_t p, const std::string varname, const unsigned degree);
    IMPLEMENT_TYPEID(URATPSERIESFLINT)
    virtual int compare(const Basic &o) const;
    virtual std::size_t __hash__() const;
    static RCP<const URatPSeriesFlint> series(const RCP<const Basic> &t, const std::string &x, unsigned int prec);
    static flint::fmpzxx convert(const Integer &x);
    static flint::fmpqxx convert(const mpq_class &x);
    static fp_t var(const std::string &s);
    static flint::fmpqxx convert(const Rational &x);
    static flint::fmpqxx convert(const Number &x);
    static inline fp_t mul(const fp_t &s, const fp_t &r, unsigned prec) {
        return fp_t(flint::mullow(s, r, prec));
    }
    static fp_t pow(const fp_t &s, int n, unsigned prec);
    static unsigned ldegree(const fp_t &s);
    static inline flint::fmpqxx find_cf(const fp_t &s, const fp_t &var, unsigned deg) {
        return flint::fmpqxx((s.get_coeff(deg)));
    }
    static flint::fmpqxx root(flint::fmpqxx &c, unsigned n);
    static fp_t diff(const fp_t &s, const fp_t &var);
    static fp_t integrate(const fp_t &s, const fp_t &var);
    static fp_t subs(const fp_t &s, const fp_t &var, const fp_t &r, unsigned prec);

    static inline fp_t series_invert(const fp_t &s, const fp_t& var, unsigned int   prec) {
        return fp_t(flint::inv_series(s, prec));
    }
    static inline fp_t series_reverse(const fp_t &s, const fp_t& var, unsigned int   prec) {
        return fp_t(flint::revert_series(s, prec));
    }
    static inline fp_t series_log(const fp_t &s, const fp_t& var, unsigned int   prec) {
        return fp_t(log_series(s, prec));
    }
    static inline fp_t series_exp(const fp_t &s, const fp_t& var, unsigned int   prec) {
        return fp_t(exp_series(s, prec));
    }
};

} // SymEngine

#endif // HAVE_SYMENGINE_FLINT

#endif //SYMENGINE_SERIES_FLINT_H
