/*
 *       /\        Guillaume Delay 2018,2019
 *      /__\       guillaume.delay@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    This is ProtoN, a library for fast Prototyping of
 *  /_\/_\/_\/_\   Numerical methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code or parts of it for scientific publications, you
 * are required to cite it as following:
 *
 * Implementation of Discontinuous Skeletal methods on arbitrary-dimensional,
 * polytopal meshes using generic programming.
 * M. Cicuttin, D. A. Di Pietro, A. Ern.
 * Journal of Computational and Applied Mathematics.
 * DOI: 10.1016/j.cam.2017.09.017
 */


template<typename T>
struct level_set {
   public:
    virtual T operator()(const disk::point<T,2>& pt) const {
    }

    virtual Eigen::Matrix<T,2,1> gradient(const disk::point<T,2>& pt) const {
    }

    Eigen::Matrix<T,2,1> normal(const disk::point<T,2>& pt) const {
        Eigen::Matrix<T,2,1> ret;
        ret = gradient(pt);
        return ret/ret.norm();
    }
};

template<typename T>
struct circle_level_set: public level_set<T>
{
    T radius, alpha, beta;

    circle_level_set(T r, T a, T b) : radius(r), alpha(a), beta(b) {
    }

    T operator()(const disk::point<T,2>& pt) const {
        auto x = pt.x();
        auto y = pt.y();
        return (x-alpha)*(x-alpha) + (y-beta)*(y-beta) - radius*radius;
    }

    Eigen::Matrix<T,2,1> gradient(const disk::point<T,2>& pt) const {
        Eigen::Matrix<T,2,1> ret;
        ret(0) = 2*pt.x() - 2*alpha;
        ret(1) = 2*pt.y() - 2*beta;
        return ret;
    }
};

template<typename T>
struct line_level_set: public level_set<T>
{
    T cut_y;

    line_level_set(T cy)
        : cut_y(cy)
    {}

    T operator()(const disk::point<T,2>& pt) const
    {
        auto x = pt.x();
        auto y = pt.y();

        return y - cut_y;
    }

    Eigen::Matrix<T,2,1> gradient(const disk::point<T,2>& pt) const {
        Eigen::Matrix<T,2,1> ret;
        ret(0) = 0;
        ret(1) = 1;
        return ret;
    }
};



template<typename T>
struct square_level_set: public level_set<T>
{
    public:
    T y_top, y_bot, x_left, x_right;

    square_level_set(T yt, T yb, T xl, T xr)
        : y_top(yt), y_bot(yb), x_left(xl), x_right(xr)
    {}

    T operator()(const disk::point<T,2>& pt) const
    {
        auto x = pt.x();
        auto y = pt.y();

        T in = 1;
        if(x > x_left && x < x_right && y > y_bot && y < y_top)
            in = 1;
        else
            in = -1;

        T dist_x = std::min( abs(x-x_left), abs(x-x_right));
        T dist_y = std::min( abs(y-y_bot), abs(y-y_top));

        
        return - in * std::min(dist_x , dist_y);
    }

    Eigen::Matrix<T,2,1> gradient(const disk::point<T,2>& pt) const
    {
        Eigen::Matrix<T,2,1> ret;
        

        auto x = pt.x();
        auto y = pt.y();

        T dist = abs(x - x_left);
        ret(0) = -1;
        ret(1) = 0;
        
        if(abs(x - x_right) < dist )
        {
            dist = abs(x - x_right);
            ret(0) = 1;
            ret(1) = 0;
        }
        if(abs(y - y_bot) < dist )
        {
            dist = abs(y - y_bot);
            ret(0) = 0;
            ret(1) = -1;
        }
        if(abs(y - y_top) < dist)
        {
            ret(0) = 0;
            ret(1) = 1;
        }
        
        return ret;
    }
};



template<typename T>
struct flower_level_set: public level_set<T>
{
    T radius, alpha, beta, a;
    size_t N;

    flower_level_set(T r, T al, T b, size_t N_, T a_)
        : radius(r), alpha(al), beta(b), N(N_), a(a_)
    {}

    T operator()(const disk::point<T,2>& pt) const
    {
        auto x = pt.x();
        auto y = pt.y();

        T theta;
        if(x == alpha && y < beta)
            theta = - M_PI / 2.0;
        else if(x == alpha && y >= beta)
            theta = M_PI / 2.0;
        else
            theta = atan((y-beta)/(x-alpha));

        if(x < alpha)
            theta = theta + M_PI;

        return (x-alpha)*(x-alpha) + (y-beta)*(y-beta) - radius*radius
            - a * std::cos(N*theta);
    }

    Eigen::Matrix<T,2,1> gradient(const disk::point<T,2>& pt) const
    {
        Eigen::Matrix<T,2,1> ret;
        auto X = pt.x() - alpha;
        auto Y = pt.y() - beta;

        T theta;
        if(X == 0 && Y < 0)
            theta = - M_PI / 2.0;
        else if(X == 0 && Y >= 0)
            theta = M_PI / 2.0;
        else
            theta = atan( Y / X );

        if(pt.x() < alpha)
            theta = theta + M_PI;
        
        ret(0) = 2*X - a * N * std::sin(N * theta) * Y / (X*X + Y*Y);
        ret(1) = 2*Y + a * N * std::sin(N * theta) * X / (X*X + Y*Y);
        return ret;
    }
};
