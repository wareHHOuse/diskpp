#pragma once

#include <vector>

#include "diskpp/quadratures/quadrature_point.hpp"

/* See mkgauss.cpp in apps/utils */

namespace disk::quadrature {

namespace priv
{

struct gauss_point
{
    double point;
    double weight;
};

static gauss_point gauss_rule_1[] = {{0, 2}};

static gauss_point gauss_rule_2[] = {{-0.577350269189626, 1}, {0.577350269189626, 1}};

static gauss_point gauss_rule_3[] = {{-0.774596669241483, 0.555555555555556},
                                     {0, 0.888888888888889},
                                     {0.774596669241483, 0.555555555555556}};

static gauss_point gauss_rule_4[] = {{-0.861136311594053, 0.347854845137454},
                                     {-0.339981043584856, 0.652145154862546},
                                     {0.339981043584856, 0.652145154862546},
                                     {0.861136311594052, 0.347854845137453}};

static gauss_point gauss_rule_5[] = {{-0.906179845938664, 0.236926885056189},
                                     {-0.538469310105683, 0.478628670499366},
                                     {0, 0.568888888888889},
                                     {0.538469310105683, 0.478628670499367},
                                     {0.906179845938664, 0.236926885056189}};

static gauss_point gauss_rule_6[] = {{-0.932469514203152, 0.17132449237917},
                                     {-0.661209386466264, 0.360761573048138},
                                     {-0.238619186083197, 0.467913934572691},
                                     {0.238619186083197, 0.467913934572691},
                                     {0.661209386466264, 0.360761573048138},
                                     {0.932469514203152, 0.171324492379171}};

static gauss_point gauss_rule_7[] = {{-0.949107912342759, 0.12948496616887},
                                     {-0.741531185599394, 0.279705391489276},
                                     {-0.405845151377397, 0.381830050505119},
                                     {0, 0.417959183673469},
                                     {0.405845151377397, 0.381830050505119},
                                     {0.741531185599394, 0.279705391489276},
                                     {0.949107912342758, 0.12948496616887}};

static gauss_point gauss_rule_8[] = {{-0.960289856497536, 0.101228536290376},
                                     {-0.796666477413627, 0.222381034453375},
                                     {-0.525532409916329, 0.313706645877887},
                                     {-0.18343464249565, 0.362683783378362},
                                     {0.18343464249565, 0.362683783378362},
                                     {0.525532409916329, 0.313706645877888},
                                     {0.796666477413627, 0.222381034453374},
                                     {0.960289856497536, 0.101228536290377}};

static gauss_point gauss_rule_9[] = {{-0.968160239507626, 0.0812743883615746},
                                     {-0.836031107326636, 0.180648160694858},
                                     {-0.61337143270059, 0.260610696402936},
                                     {-0.324253423403809, 0.312347077040003},
                                     {0, 0.33023935500126},
                                     {0.324253423403809, 0.312347077040003},
                                     {0.613371432700591, 0.260610696402936},
                                     {0.836031107326637, 0.180648160694857},
                                     {0.968160239507627, 0.0812743883615752}};

static gauss_point gauss_rule_10[] = {{-0.973906528517171, 0.0666713443086882},
                                      {-0.865063366688984, 0.14945134915058},
                                      {-0.679409568299024, 0.219086362515982},
                                      {-0.433395394129247, 0.269266719309996},
                                      {-0.148874338981631, 0.295524224714753},
                                      {0.148874338981631, 0.295524224714753},
                                      {0.433395394129247, 0.269266719309997},
                                      {0.679409568299025, 0.219086362515982},
                                      {0.865063366688985, 0.149451349150582},
                                      {0.973906528517172, 0.0666713443086878}};

static gauss_point gauss_rule_11[] = {{-0.978228658146058, 0.0556685671161741},
                                      {-0.887062599768096, 0.125580369464904},
                                      {-0.73015200557405, 0.186290210927735},
                                      {-0.519096129206812, 0.233193764591991},
                                      {-0.269543155952345, 0.262804544510247},
                                      {0, 0.272925086777901},
                                      {0.269543155952345, 0.262804544510247},
                                      {0.519096129206813, 0.233193764591992},
                                      {0.73015200557405, 0.186290210927734},
                                      {0.887062599768097, 0.125580369464905},
                                      {0.978228658146058, 0.0556685671161723}};

static gauss_point gauss_rule_12[] = {{-0.981560634246719, 0.0471753363865106},
                                      {-0.904117256370475, 0.106939325995319},
                                      {-0.769902674194305, 0.160078328543347},
                                      {-0.587317954286618, 0.203167426723066},
                                      {-0.36783149899818, 0.233492536538355},
                                      {-0.125233408511469, 0.249147045813402},
                                      {0.125233408511469, 0.249147045813403},
                                      {0.36783149899818, 0.233492536538355},
                                      {0.587317954286617, 0.203167426723066},
                                      {0.769902674194304, 0.160078328543346},
                                      {0.904117256370473, 0.106939325995317},
                                      {0.981560634246718, 0.047175336386512}};

static gauss_point gauss_rule_13[] = {{-0.984183054718589, 0.0404840047653166},
                                      {-0.917598399222978, 0.0921214998377284},
                                      {-0.80157809073331, 0.138873510219787},
                                      {-0.64234933944034, 0.178145980761946},
                                      {-0.448492751036447, 0.207816047536889},
                                      {-0.230458315955135, 0.226283180262897},
                                      {0, 0.232551553230875},
                                      {0.230458315955135, 0.226283180262897},
                                      {0.448492751036447, 0.207816047536888},
                                      {0.64234933944034, 0.178145980761946},
                                      {0.80157809073331, 0.138873510219787},
                                      {0.917598399222979, 0.0921214998377297},
                                      {0.984183054718589, 0.0404840047653156}};

static gauss_point gauss_rule_14[] = {{-0.986283808696813, 0.0351194603317527},
                                      {-0.928434883663573, 0.0801580871597601},
                                      {-0.827201315069765, 0.121518570687902},
                                      {-0.687292904811685, 0.157203167158194},
                                      {-0.515248636358154, 0.185538397477938},
                                      {-0.31911236892789, 0.205198463721296},
                                      {-0.108054948707344, 0.215263853463158},
                                      {0.108054948707344, 0.215263853463158},
                                      {0.31911236892789, 0.205198463721296},
                                      {0.515248636358154, 0.185538397477938},
                                      {0.687292904811686, 0.157203167158195},
                                      {0.827201315069764, 0.121518570687902},
                                      {0.928434883663573, 0.0801580871597608},
                                      {0.986283808696812, 0.035119460331751}};

static gauss_point gauss_rule_15[] = {{-0.987992518020485, 0.0307532419961169},
                                      {-0.937273392400705, 0.0703660474881079},
                                      {-0.848206583410427, 0.107159220467172},
                                      {-0.72441773136017, 0.139570677926155},
                                      {-0.570972172608539, 0.166269205816994},
                                      {-0.394151347077564, 0.186161000015562},
                                      {-0.201194093997435, 0.198431485327112},
                                      {0, 0.202578241925561},
                                      {0.201194093997435, 0.198431485327112},
                                      {0.394151347077563, 0.186161000015562},
                                      {0.570972172608538, 0.166269205816994},
                                      {0.72441773136017, 0.139570677926155},
                                      {0.848206583410426, 0.107159220467171},
                                      {0.937273392400705, 0.0703660474881086},
                                      {0.987992518020484, 0.0307532419961168}};

static gauss_point gauss_rule_16[] = {{-0.98940093499165, 0.0271524594117561},
                                      {-0.944575023073232, 0.0622535239386456},
                                      {-0.865631202387832, 0.0951585116824928},
                                      {-0.755404408355003, 0.124628971255534},
                                      {-0.617876244402644, 0.149595988816577},
                                      {-0.458016777657227, 0.169156519395002},
                                      {-0.281603550779259, 0.182603415044923},
                                      {-0.0950125098376372, 0.189450610455068},
                                      {0.0950125098376373, 0.189450610455068},
                                      {0.281603550779259, 0.182603415044923},
                                      {0.458016777657227, 0.169156519395003},
                                      {0.617876244402644, 0.149595988816576},
                                      {0.755404408355004, 0.124628971255536},
                                      {0.865631202387832, 0.0951585116824934},
                                      {0.944575023073233, 0.0622535239386473},
                                      {0.98940093499165, 0.0271524594117542}};

static gauss_point gauss_rule_17[] = {{-0.990575475314417, 0.0241483028685482},
                                      {-0.950675521768768, 0.0554595293739856},
                                      {-0.880239153726986, 0.0850361483171798},
                                      {-0.781514003896801, 0.111883847193404},
                                      {-0.657671159216691, 0.135136368468526},
                                      {-0.512690537086477, 0.15404576107681},
                                      {-0.351231763453876, 0.16800410215645},
                                      {-0.178484181495848, 0.176562705366993},
                                      {0, 0.179446470356207},
                                      {0.178484181495848, 0.176562705366992},
                                      {0.351231763453876, 0.16800410215645},
                                      {0.512690537086476, 0.154045761076811},
                                      {0.65767115921669, 0.135136368468525},
                                      {0.781514003896801, 0.111883847193405},
                                      {0.880239153726984, 0.085036148317177},
                                      {0.950675521768767, 0.0554595293739875},
                                      {0.990575475314416, 0.0241483028685473}};

static gauss_point gauss_rule_18[] = {{-0.991565168420931, 0.021616013526483},
                                      {-0.955823949571398, 0.0497145488949712},
                                      {-0.892602466497556, 0.0764257302548879},
                                      {-0.803704958972523, 0.100942044106287},
                                      {-0.691687043060353, 0.122555206711478},
                                      {-0.559770831073947, 0.140642914670651},
                                      {-0.411751161462843, 0.154684675126265},
                                      {-0.251886225691506, 0.164276483745832},
                                      {-0.0847750130417353, 0.169142382963144},
                                      {0.0847750130417352, 0.169142382963143},
                                      {0.251886225691506, 0.164276483745833},
                                      {0.411751161462843, 0.154684675126266},
                                      {0.559770831073947, 0.14064291467065},
                                      {0.691687043060353, 0.122555206711478},
                                      {0.803704958972522, 0.100942044106286},
                                      {0.892602466497556, 0.0764257302548886},
                                      {0.955823949571398, 0.049714548894969},
                                      {0.991565168420932, 0.0216160135264855}};

static gauss_point gauss_rule_19[] = {{-0.992406843843583, 0.0194617882297263},
                                      {-0.96020815213483, 0.0448142267656976},
                                      {-0.903155903614819, 0.0690445427376437},
                                      {-0.822714656537143, 0.0914900216224497},
                                      {-0.72096617733523, 0.111566645547334},
                                      {-0.60054530466168, 0.128753962539336},
                                      {-0.464570741375961, 0.142606702173606},
                                      {-0.31656409996363, 0.15276604206586},
                                      {-0.160358645640225, 0.158968843393954},
                                      {0, 0.161054449848783},
                                      {0.160358645640225, 0.158968843393954},
                                      {0.31656409996363, 0.15276604206586},
                                      {0.46457074137596, 0.142606702173606},
                                      {0.600545304661681, 0.128753962539336},
                                      {0.720966177335229, 0.111566645547334},
                                      {0.822714656537143, 0.0914900216224515},
                                      {0.903155903614818, 0.0690445427376399},
                                      {0.960208152134831, 0.0448142267656984},
                                      {0.992406843843585, 0.0194617882297286}};

static gauss_point gauss_rule_20[] = {
  {-0.993128599185095, 0.0176140071391534}, {-0.963971927277914, 0.0406014298003863},
  {-0.912234428251325, 0.0626720483341077}, {-0.839116971822219, 0.0832767415767044},
  {-0.746331906460151, 0.101930119817241},  {-0.636053680726515, 0.118194531961519},
  {-0.510867001950827, 0.131688638449177},  {-0.37370608871542, 0.142096109318382},
  {-0.227785851141645, 0.149172986472604},  {-0.0765265211334973, 0.152753387130726},
  {0.0765265211334971, 0.152753387130725},  {0.227785851141645, 0.149172986472603},
  {0.37370608871542, 0.142096109318383},    {0.510867001950828, 0.131688638449177},
  {0.636053680726517, 0.118194531961519},   {0.746331906460152, 0.10193011981724},
  {0.839116971822218, 0.0832767415767037},  {0.912234428251327, 0.0626720483341074},
  {0.963971927277915, 0.0406014298003887},  {0.993128599185096, 0.017614007139152}};

struct gauss_rule
{
    size_t       num_entries;
    gauss_point* points;
};

static struct gauss_rule gauss_rules[] = {{sizeof(gauss_rule_1) / (sizeof(gauss_point)), gauss_rule_1},
                                          {sizeof(gauss_rule_2) / (sizeof(gauss_point)), gauss_rule_2},
                                          {sizeof(gauss_rule_3) / (sizeof(gauss_point)), gauss_rule_3},
                                          {sizeof(gauss_rule_4) / (sizeof(gauss_point)), gauss_rule_4},
                                          {sizeof(gauss_rule_5) / (sizeof(gauss_point)), gauss_rule_5},
                                          {sizeof(gauss_rule_6) / (sizeof(gauss_point)), gauss_rule_6},
                                          {sizeof(gauss_rule_7) / (sizeof(gauss_point)), gauss_rule_7},
                                          {sizeof(gauss_rule_8) / (sizeof(gauss_point)), gauss_rule_8},
                                          {sizeof(gauss_rule_9) / (sizeof(gauss_point)), gauss_rule_9},
                                          {sizeof(gauss_rule_10) / (sizeof(gauss_point)), gauss_rule_10},
                                          {sizeof(gauss_rule_11) / (sizeof(gauss_point)), gauss_rule_11},
                                          {sizeof(gauss_rule_12) / (sizeof(gauss_point)), gauss_rule_12},
                                          {sizeof(gauss_rule_13) / (sizeof(gauss_point)), gauss_rule_13},
                                          {sizeof(gauss_rule_14) / (sizeof(gauss_point)), gauss_rule_14},
                                          {sizeof(gauss_rule_15) / (sizeof(gauss_point)), gauss_rule_15},
                                          {sizeof(gauss_rule_16) / (sizeof(gauss_point)), gauss_rule_16},
                                          {sizeof(gauss_rule_17) / (sizeof(gauss_point)), gauss_rule_17},
                                          {sizeof(gauss_rule_18) / (sizeof(gauss_point)), gauss_rule_18},
                                          {sizeof(gauss_rule_19) / (sizeof(gauss_point)), gauss_rule_19},
                                          {sizeof(gauss_rule_20) / (sizeof(gauss_point)), gauss_rule_20},
                                          {0, NULL}};

#define QUAD_RAW_GAUSS_LEGENDRE_MAX_ORDER 39

} // namespace priv

template<typename T>
std::vector<quadrature_point<T, 1>>
gauss_legendre(size_t degree, const T& a, const T& b)
{
    const auto num_rules = sizeof(priv::gauss_rules) / sizeof(priv::gauss_rule) - 1;
    const auto rule_num  = degree / 2;
    if (rule_num >= num_rules)
        throw std::invalid_argument("gauss_legendre: order too high");

    auto                                npts = priv::gauss_rules[rule_num].num_entries;
    std::vector<quadrature_point<T, 1>> qps;
    qps.reserve(npts);
    for (size_t i = 0; i < npts; i++) {
        const auto &qp = priv::gauss_rules[rule_num].points[i];
        auto tr_qp = 0.5*(a+b) + 0.5*(b-a)*qp.point;
        auto tr_qw = 0.5 * (b-a) * qp.weight;
        qps.push_back({tr_qp, tr_qw});
    }

    return qps;
}

template<typename T, size_t DIM>
std::vector<quadrature_point<T, DIM>>
gauss_legendre(size_t degree, const point<T, DIM>& a, const point<T, DIM>& b)
{
    const auto num_rules = sizeof(priv::gauss_rules) / sizeof(priv::gauss_rule) - 1;
    const auto rule_num  = degree / 2;
    if (rule_num >= num_rules)
        throw std::invalid_argument("gauss_legendre: order too high");

    const auto&                           qrule = priv::gauss_rules[rule_num];
    auto                                  npts  = qrule.num_entries;
    std::vector<quadrature_point<T, DIM>> qps;
    qps.reserve(npts);
    for (size_t i = 0; i < npts; i++)
    {
        const auto& qp    = qrule.points[i];
        auto        tr_qp = 0.5 * (a + b) + 0.5 * (b - a) * qp.point;
        auto        tr_qw = 0.5 * distance(a, b) * qp.weight;
        qps.push_back({tr_qp, tr_qw});
    }

    return qps;
}

template<typename T>
std::vector<quadrature_point<T, 2>>
tensorized_gauss_legendre(size_t             degree,
                          const point<T, 2>& p0,
                          const point<T, 2>& p1,
                          const point<T, 2>& p2,
                          const point<T, 2>& p3)
{
    const auto num_rules = sizeof(priv::gauss_rules) / sizeof(priv::gauss_rule) - 1;
    const auto rule_num  = degree / 2;
    if (rule_num >= num_rules)
        throw std::invalid_argument("tensorized_gauss_legendre: order too high");

    const auto& qrule = priv::gauss_rules[rule_num];

    auto                                npts = qrule.num_entries;
    std::vector<quadrature_point<T, 2>> qps;
    qps.reserve(npts * npts); /* Collect the points in the [-1, 1]^2 square */
    for (size_t i = 0; i < npts; i++)
    {
        const auto& qpy = qrule.points[i];
        for (size_t j = 0; j < npts; j++)
        {
            const auto& qpx = qrule.points[j];
            point<T, 2> qp(qpx.point, qpy.point);
            auto        qw = qpx.weight * qpy.weight;
            qps.push_back({qp, qw});
        }
    }

    /* points in pts must be in _counterclockwise_ order */
    auto X = [&](T xi, T eta) -> T
    {
        return 0.25 * p0.x() * (1 - xi) * (1 - eta) + 0.25 * p1.x() * (1 + xi) * (1 - eta) +
               0.25 * p2.x() * (1 + xi) * (1 + eta) + 0.25 * p3.x() * (1 - xi) * (1 + eta);
    };

    auto Y = [&](T xi, T eta) -> T
    {
        return 0.25 * p0.y() * (1 - xi) * (1 - eta) + 0.25 * p1.y() * (1 + xi) * (1 - eta) +
               0.25 * p2.y() * (1 + xi) * (1 + eta) + 0.25 * p3.y() * (1 - xi) * (1 + eta);
    };

    auto J = [&](T xi, T eta) -> T
    {
        auto j11 = 0.25 * ((p1.x() - p0.x()) * (1 - eta) + (p2.x() - p3.x()) * (1 + eta));

        auto j12 = 0.25 * ((p1.y() - p0.y()) * (1 - eta) + (p2.y() - p3.y()) * (1 + eta));

        auto j21 = 0.25 * ((p3.x() - p0.x()) * (1 - xi) + (p2.x() - p1.x()) * (1 + xi));

        auto j22 = 0.25 * ((p3.y() - p0.y()) * (1 - xi) + (p2.y() - p1.y()) * (1 + xi));

        return std::abs(j11 * j22 - j12 * j21);
    };

    /* do ref->phys transform in place */
    for (auto& qp : qps)
    {
        const T xi  = qp.point().x();
        const T eta = qp.point().y();

        point<T, 2> phys_point(X(xi, eta), Y(xi, eta));
        const T     phys_weight = qp.weight() * J(xi, eta);
        qp                      = disk::make_qp(phys_point, phys_weight);
    }

    return qps;
}

template<typename T>
std::vector<quadrature_point<T, 2>>
tensorized_gauss_legendre(const size_t degree, const std::array<point<T, 2>, 4>& pts)
{
    return tensorized_gauss_legendre(degree, pts[0], pts[1], pts[2], pts[3]);
}

template<typename T>
std::vector<quadrature_point<T, 3>>
tensorized_gauss_legendre(size_t             degree,
                          const point<T, 3>& p0,
                          const point<T, 3>& p1,
                          const point<T, 3>& p2,
                          const point<T, 3>& p3)
{
    const auto num_rules = sizeof(priv::gauss_rules) / sizeof(priv::gauss_rule) - 1;
    const auto rule_num  = degree / 2;
    if (rule_num >= num_rules)
        throw std::invalid_argument("tensorized_gauss_legendre: order too high");

    const auto& qrule = priv::gauss_rules[rule_num];

    auto                                npts = qrule.num_entries;
    std::vector<quadrature_point<T, 2>> qps_2d;
    qps_2d.reserve(npts * npts); /* Collect the points in the [-1, 1]^2 square */
    for (size_t i = 0; i < npts; i++)
    {
        const auto& qpy = qrule.points[i];
        for (size_t j = 0; j < npts; j++)
        {
            const auto& qpx = qrule.points[j];
            point<T, 2> qp(qpx.point, qpy.point);
            auto        qw = qpx.weight * qpy.weight;
            qps_2d.push_back({qp, qw});
        }
    }

    const auto v0    = (p1 - p0) / T(2);
    const auto v1    = (p3 - p0) / T(2);
    const T    meas4 = v0.norm()*v1.norm();

    const auto bar = (p0 + p1 + p2 + p3) / T(4);

    /* do ref->phys transform in place */
    std::vector<quadrature_point<T, 3>> qps;
    qps.reserve(qps_2d.size());
    for (auto& qp : qps_2d)
    {
        const auto qpc    = qp.point();
        const auto point  = v0 * qpc.x() + v1 * qpc.y() + bar;
        const auto weight = qp.weight() * meas4;
        qps.push_back(disk::make_qp(point, weight));
    }

    return qps;
}

template<typename T>
std::vector<disk::quadrature_point<T, 3>>
tensorized_gauss_legendre(const size_t degree, const std::array<point<T, 3>, 8>& pts)
{
    const auto num_rules = sizeof(priv::gauss_rules) / sizeof(priv::gauss_rule) - 1;
    const auto rule_num  = degree / 2;
    if (rule_num >= num_rules)
        throw std::invalid_argument("gauss_legendre: order too high");

    const auto&                           qrule = priv::gauss_rules[rule_num];
    auto                                  npts  = qrule.num_entries;
    std::vector<quadrature_point<T, 1>> qps;
    qps.reserve(npts);
    for (size_t i = 0; i < npts; i++)
    {
        const auto& qp    = qrule.points[i];
        qps.push_back({qp.point, qp.weight});
    }

    // std::cout << "deg:" << degree << " nb : " << 3 * qps.size() << std::endl;

    std::vector<disk::quadrature_point<T, 3>> ret;
    ret.reserve(3 * qps.size());

    auto P = [&pts](T xi, T eta, T zeta) -> T
    {
        const T val = pts[0].x() * (1 - xi) * (1 - eta) * (1 - zeta) + pts[1].x() * (1 + xi) * (1 - eta) * (1 - zeta) +
                      pts[2].x() * (1 + xi) * (1 + eta) * (1 - zeta) + pts[3].x() * (1 - xi) * (1 + eta) * (1 - zeta) +
                      pts[4].x() * (1 - xi) * (1 - eta) * (1 + zeta) + pts[5].x() * (1 + xi) * (1 - eta) * (1 + zeta) +
                      pts[6].x() * (1 + xi) * (1 + eta) * (1 + zeta) + pts[7].x() * (1 - xi) * (1 + eta) * (1 + zeta);
        return T(0.125) * val;
    };

    auto Q = [&pts](T xi, T eta, T zeta) -> T
    {
        const T val = pts[0].y() * (1 - xi) * (1 - eta) * (1 - zeta) + pts[1].y() * (1 + xi) * (1 - eta) * (1 - zeta) +
                      pts[2].y() * (1 + xi) * (1 + eta) * (1 - zeta) + pts[3].y() * (1 - xi) * (1 + eta) * (1 - zeta) +
                      pts[4].y() * (1 - xi) * (1 - eta) * (1 + zeta) + pts[5].y() * (1 + xi) * (1 - eta) * (1 + zeta) +
                      pts[6].y() * (1 + xi) * (1 + eta) * (1 + zeta) + pts[7].y() * (1 - xi) * (1 + eta) * (1 + zeta);
        return T(0.125) * val;
    };

    auto R = [&pts](T xi, T eta, T zeta) -> T
    {
        const T val = pts[0].z() * (1 - xi) * (1 - eta) * (1 - zeta) + pts[1].z() * (1 + xi) * (1 - eta) * (1 - zeta) +
                      pts[2].z() * (1 + xi) * (1 + eta) * (1 - zeta) + pts[3].z() * (1 - xi) * (1 + eta) * (1 - zeta) +
                      pts[4].z() * (1 - xi) * (1 - eta) * (1 + zeta) + pts[5].z() * (1 + xi) * (1 - eta) * (1 + zeta) +
                      pts[6].z() * (1 + xi) * (1 + eta) * (1 + zeta) + pts[7].z() * (1 - xi) * (1 + eta) * (1 + zeta);
        return T(0.125) * val;
    };

    auto J = [&pts](T xi, T eta, T zeta) -> T
    {
        static_matrix<T, 3, 3> Jac = static_matrix<T, 3, 3>::Zero();
        Jac(0, 0)                  = -pts[0].x() * (1 - eta) * (1 - zeta) + pts[1].x() * (1 - eta) * (1 - zeta) +
                    pts[2].x() * (1 + eta) * (1 - zeta) - pts[3].x() * (1 + eta) * (1 - zeta) -
                    pts[4].x() * (1 - eta) * (1 + zeta) + pts[5].x() * (1 - eta) * (1 + zeta) +
                    pts[6].x() * (1 + eta) * (1 + zeta) - pts[7].x() * (1 + eta) * (1 + zeta);

        Jac(0, 1) = -pts[0].y() * (1 - eta) * (1 - zeta) + pts[1].y() * (1 - eta) * (1 - zeta) +
                    pts[2].y() * (1 + eta) * (1 - zeta) - pts[3].y() * (1 + eta) * (1 - zeta) -
                    pts[4].y() * (1 - eta) * (1 + zeta) + pts[5].y() * (1 - eta) * (1 + zeta) +
                    pts[6].y() * (1 + eta) * (1 + zeta) - pts[7].y() * (1 + eta) * (1 + zeta);

        Jac(0, 2) = -pts[0].z() * (1 - eta) * (1 - zeta) + pts[1].z() * (1 - eta) * (1 - zeta) +
                    pts[2].z() * (1 + eta) * (1 - zeta) - pts[3].z() * (1 + eta) * (1 - zeta) -
                    pts[4].z() * (1 - eta) * (1 + zeta) + pts[5].z() * (1 - eta) * (1 + zeta) +
                    pts[6].z() * (1 + eta) * (1 + zeta) - pts[7].z() * (1 + eta) * (1 + zeta);

        Jac(1, 0) = -pts[0].x() * (1 - xi) * (1 - zeta) - pts[1].x() * (1 + xi) * (1 - zeta) +
                    pts[2].x() * (1 + xi) * (1 - zeta) + pts[3].x() * (1 - xi) * (1 - zeta) -
                    pts[4].x() * (1 - xi) * (1 + zeta) - pts[5].x() * (1 + xi) * (1 + zeta) +
                    pts[6].x() * (1 + xi) * (1 + zeta) + pts[7].x() * (1 - xi) * (1 + zeta);

        Jac(1, 1) = -pts[0].y() * (1 - xi) * (1 - zeta) - pts[1].y() * (1 + xi) * (1 - zeta) +
                    pts[2].y() * (1 + xi) * (1 - zeta) + pts[3].y() * (1 - xi) * (1 - zeta) -
                    pts[4].y() * (1 - xi) * (1 + zeta) - pts[5].y() * (1 + xi) * (1 + zeta) +
                    pts[6].y() * (1 + xi) * (1 + zeta) + pts[7].y() * (1 - xi) * (1 + zeta);

        Jac(1, 2) = -pts[0].z() * (1 - xi) * (1 - zeta) - pts[1].z() * (1 + xi) * (1 - zeta) +
                    pts[2].z() * (1 + xi) * (1 - zeta) + pts[3].z() * (1 - xi) * (1 - zeta) -
                    pts[4].z() * (1 - xi) * (1 + zeta) - pts[5].z() * (1 + xi) * (1 + zeta) +
                    pts[6].z() * (1 + xi) * (1 + zeta) + pts[7].z() * (1 - xi) * (1 + zeta);

        Jac(2, 0) = -pts[0].x() * (1 - xi) * (1 - eta) - pts[1].x() * (1 + xi) * (1 - eta) -
                    pts[2].x() * (1 + xi) * (1 + eta) - pts[3].x() * (1 - xi) * (1 + eta) +
                    pts[4].x() * (1 - xi) * (1 - eta) + pts[5].x() * (1 + xi) * (1 - eta) +
                    pts[6].x() * (1 + xi) * (1 + eta) + pts[7].x() * (1 - xi) * (1 + eta);

        Jac(2, 1) = -pts[0].y() * (1 - xi) * (1 - eta) - pts[1].y() * (1 + xi) * (1 - eta) -
                    pts[2].y() * (1 + xi) * (1 + eta) - pts[3].y() * (1 - xi) * (1 + eta) +
                    pts[4].y() * (1 - xi) * (1 - eta) + pts[5].y() * (1 + xi) * (1 - eta) +
                    pts[6].y() * (1 + xi) * (1 + eta) + pts[7].y() * (1 - xi) * (1 + eta);

        Jac(2, 2) = -pts[0].z() * (1 - xi) * (1 - eta) - pts[1].z() * (1 + xi) * (1 - eta) -
                    pts[2].z() * (1 + xi) * (1 + eta) - pts[3].z() * (1 - xi) * (1 + eta) +
                    pts[4].z() * (1 - xi) * (1 - eta) + pts[5].z() * (1 + xi) * (1 - eta) +
                    pts[6].z() * (1 + xi) * (1 + eta) + pts[7].z() * (1 - xi) * (1 + eta);

        Jac *= T(0.125);
        return std::abs(Jac.determinant());
    };

    for (auto& qp_x : qps)
    {
        const auto xi = qp_x.point().x();
        const auto w_x   = qp_x.weight();
        for (auto& qp_y : qps)
        {
            const auto eta   = qp_y.point().x();
            const auto w_y   = qp_y.weight();
            for (auto& qp_z : qps)
            {
                const auto zeta  = qp_z.point().x();
                const auto w_z   = qp_z.weight();

                const auto px = P(xi, eta, zeta);
                const auto py = Q(xi, eta, zeta);
                const auto pz = R(xi, eta, zeta);

                // std::cout << xi << " " << px << std::endl;
                // std::cout << eta << " " << py << std::endl;
                // std::cout << zeta << " " << pz << std::endl;
                // std::cout << J(xi, eta, zeta) << std::endl;

                const auto w = w_x * w_y * w_z * J(xi, eta, zeta);

                ret.push_back(disk::make_qp(point<T, 3>({px, py, pz}), w));
            }
        }
    }

    return ret;
}

} // namespace disk::quadrature