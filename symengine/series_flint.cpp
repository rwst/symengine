
#include <symengine/series_flint.h>
#include <symengine/series_visitor.h>

namespace SymEngine {

URatPSeriesFlint::URatPSeriesFlint(fp_t p, const std::string varname, const unsigned degree)
        : SeriesBase(std::move(p), varname, degree) {

}
RCP<const URatPSeriesFlint> URatPSeriesFlint::series(const RCP<const Basic> &t, const std::string &x,
                                                         unsigned int prec) {
    SeriesVisitor<fp_t, flint::fmpqxx, URatPSeriesFlint, false> visitor(fp_t(x), x, prec);
    return visitor.series(t);
}

std::size_t URatPSeriesFlint::__hash__() const {
    std::size_t seed = URATPSERIESPIRANHA;
    hash_combine(seed, p_.hash());
    hash_combine(seed, var_);
    hash_combine(seed, degree_);
    return seed;
}

int URatPSeriesFlint::compare(const Basic &o) const {
    SYMENGINE_ASSERT(is_a<URatPSeriesFlint>(o))
    const URatPSeriesFlint &s = static_cast<const URatPSeriesFlint &>(o);
    if (var_ != s.var_)
        return (var_ < s.var_) ? -1 : 1;
    if (degree_ != s.degree_)
        return (degree_ < s.degree_) ? -1 : 1;
    if (p_ == s.p_)
        return 0;
    return (p_.hash() < s.p_.hash()) ? -1 : 1;
}

fmpzxx URatPSeriesFlint::convert(const Integer &x) {
    return fmpzxx(x.as_mpz().get_mpz_t());
}

fmpqxx URatPSeriesFlint::convert(const mpq_class &x) {
    fmpzxx i1(x.get_num_mpz_t());
    fmpzxx i2(x.get_den_mpz_t());
    fmpqxx r(i1);
    r /= i2;
    return r;
}

fp_t URatPSeriesFlint::var(const std::string &s) {
    return fp_t(s);
}

fmpqxx URatPSeriesFlint::convert(const Rational &x) {
    return convert(x.as_mpq());
}

fmpqxx URatPSeriesFlint::convert(const Number &x) {
    throw std::runtime_error("Not Implemented");
}

fp_t URatPSeriesFlint::pow(const fp_t &s, int n, unsigned prec) {
    if (n > 0)
        return flint::pow(s, n);
    else if (n < 0)
        return flint::pow(flint::inv_series(s, prec), -n);
    return fp_t(1);
}

unsigned URatPSeriesFlint::ldegree(const fp_t &s) {
    long i = 0;
    while (i <= s.degree())
        if (not s.get_coeff(i++).is_zero())
            return i-1;
    return 0;
}

fmpqxx URatPSeriesFlint::root(fmpqxx &c, unsigned n) {
    mpq_t g;
    mpq_init(g);
    fmpq_get_mpq(g, c);
    mpq_class cl_rat(g);
    mpq_clean();
    mpq_class cl_root;
    bool res;
    if (mpz_cmp_ui(cl_rat.get_den_mpz_t(), 1) == 0) {
        // integer constant
        res = mpz_root(cl_root.get_num_mpz_t(), cl_rat.get_num_mpz_t(), n) != 0;
    }
    else {
        cl_root.canonicalize();
        RCP<const Rational> cterm = make_rcp<const Rational>(cl_root);
        RCP<const Number> cout;
        res = cterm->nth_root(outArg(cout), n);
        cl_root = static_cast<const Rational&>(*cout).i;
    }
    if (not res)
        throw std::runtime_error("constant term is not an nth power");
    return convert(cl_root);
}

fp_t URatPSeriesFlint::diff(const fp_t &s, const fp_t &var) {
    return flint::derivative(s);
}

fp_t URatPSeriesFlint::integrate(const fp_t &s, const fp_t &var) {
    return flint::integral(s);
}

fp_t URatPSeriesFlint::subs(const fp_t &s, const fp_t &var, const fp_t &r, unsigned prec) {
    return s(r);
}

}
}
