#include <fstream>
#include "diskpp/quadratures/quadratures.hpp"
#include "diskpp/loaders/loader.hpp"
#include "diskpp/output/silo.hpp"


template<typename T, size_t DIM>
T scaling_factor(const std::vector<disk::point<T,DIM>>& points)
{
    return 0.0;
}

template<typename T>
T scaling_factor(const std::vector<disk::point<T,2>>& points)
{
    T min_x = points[0].x();
    T min_y = points[0].y();
    T max_x = points[0].x();
    T max_y = points[0].y();

    for (auto& pt : points)
    {
        min_x = std::min(min_x, pt.x());
        min_y = std::min(min_y, pt.y());
        max_x = std::max(max_x, pt.x());
        max_y = std::max(max_y, pt.y());
    }

    return 0.1*std::max( max_x-min_x, max_y-min_y );
}

template<typename T, size_t DIM>
disk::point<T,DIM>
barycenter(const std::vector<disk::point<T,DIM>>& points)
{
    disk::point<T,DIM> ret = std::accumulate(points.begin(), points.end(), disk::point<T,DIM>{});
    return ret * 1./T(points.size());
}

template<typename T, size_t DIM>
T diameter(const std::vector<disk::point<T,DIM>>& points)
{
    T diam = 0.0;
    for (size_t i = 0; i < points.size(); i++)
        for (size_t j = i; j < points.size(); j++)
            diam = std::max( diam, distance(points[i], points[j]) );
    
    return diam;
}



class postscript_plotter
{
    static constexpr double paper_width = 210.;
    static constexpr double paper_height = 297.;
    static constexpr double box_width = 50.4;
    static constexpr double box_height = 25.4;
    static constexpr double margin_x = 5;
    static constexpr double margin_y = 5;
    static constexpr int rows = paper_height/box_height;
    static constexpr int cols = paper_width/box_width;

    int curr_index;

    std::ofstream ps_ofs;

    double mm_to_points(double mm)
    {
        return 72.*mm/25.4;
    }

    double points_to_mm(double points)
    {
        return points*25.4/72.;
    }

    

    template<typename T>
    std::pair<disk::point<T,2>, disk::point<T,2>>
    centers(int cell_index)
    {
        auto num_cells = rows*cols;
        cell_index %= num_cells;

        auto row = cell_index / cols;
        auto col = cell_index % cols;

        double dx = box_width/4.0; 
        double center1_x = margin_x + col*box_width + dx;
        double center2_x = margin_x + col*box_width + 3*dx;

        double top = margin_y + rows*box_height;
        double dy = box_height/2.0;
        double center_y = top - row*box_height - dy;

        disk::point<T,2> c1(center1_x, center_y);
        disk::point<T,2> c2(center2_x, center_y);
        return std::make_pair(c1, c2);
    }

    void linewidth(double lw)
    {
        ps_ofs << lw << " setlinewidth" << std::endl;
    }

    void line(double sx, double sy, double ex, double ey)
    {
        ps_ofs << sx << " " << sy << " moveto" << std::endl;
        ps_ofs << ex << " " << ey << " lineto" << std::endl;
    }

    void moveto(double x, double y)
    {
        ps_ofs << x << " " << y << " moveto" << std::endl;
    }

    void rlineto(double x, double y)
    {
        ps_ofs << x << " " << y << " rlineto" << std::endl;
    }

    void setrgbcolor(double r, double g, double b)
    {
        ps_ofs << r << " " << g << " " << b << " setrgbcolor" << std::endl;
    }

    void closepath(void)
    {
        ps_ofs << "closepath" << std::endl;
    }

    void stroke(void)
    {
        ps_ofs << "stroke" << std::endl;
    }

    void show_page(void)
    {
        ps_ofs << "showpage" << std::endl;
    }

    template<typename T>
    void write(const T& t)
    {
        ps_ofs << "(" << t << ") show" << std::endl;
    }

    void make_grid(void)
    {
        setrgbcolor(0.,0.,0.);
        linewidth(0.2);

        for (size_t row = 0; row < rows+1; row++)
        {
            double line_sx = mm_to_points(margin_x);
            double line_sy = mm_to_points(margin_y + row*box_height);

            double line_ex = mm_to_points(margin_x + cols*box_width);
            double line_ey = mm_to_points(margin_y + row*box_height);

            line(line_sx, line_sy, line_ex, line_ey);
        }

        for (size_t col = 0; col < cols+1; col++)
        {
            double line_sx = mm_to_points(margin_x + col*box_width);
            double line_sy = mm_to_points(margin_y);

            double line_ex = mm_to_points(margin_x + col*box_width);
            double line_ey = mm_to_points(margin_y + rows*box_height);

            line(line_sx, line_sy, line_ex, line_ey);
        }

        stroke();
    }
    
public:
    postscript_plotter()
        : curr_index(0)
    {}

    postscript_plotter(const std::string& filename)
        : curr_index(0)
    {
        ps_ofs.open(filename);
        ps_ofs << "<< /PageSize [595 842] >> setpagedevice" << std::endl;
        ps_ofs << "/Courier findfont" << std::endl;
        ps_ofs << "8 scalefont setfont" << std::endl;
    }

    template<typename T>
    void add_poly(const std::vector<disk::point<T,3>>& points, const Eigen::Matrix<T,3,3>& m)
    {}

    template<typename T>
    void add_poly(const std::vector<disk::point<T,2>>& points,
        const Eigen::Matrix<T,2,2>& m)
    {
        auto tr = [&](const disk::point<T,2>& pt) {
            Eigen::Matrix<T, 2, 1> trp = m.transpose()*pt.to_vector();
            disk::point<T,2> ret;
            ret.x() = trp(0);
            ret.y() = trp(1);
            return ret;
        };

        std::vector<disk::point<T,2>> points_tr(points.size());
        std::transform(points.begin(), points.end(), points_tr.begin(), tr);

        if (curr_index == 0) {
            make_grid();
        }
        
        setrgbcolor(0.,0.,0.);
        linewidth(0.6);
        auto [c1, c2] = centers<T>(curr_index);
        
        auto sf_orig = scaling_factor(points);
        auto bar_orig = barycenter(points);
        auto p0 = c1 + (points[0] - bar_orig)/sf_orig;
        moveto( mm_to_points(p0.x()), mm_to_points(p0.y()) );
        for (size_t i = 1; i < points.size(); i++) {
            auto v = (points[i] - points[i-1])/sf_orig;
            rlineto( mm_to_points(v.x()), mm_to_points(v.y()) );
        }
        closepath();

        moveto( mm_to_points(c1.x()-8.0), mm_to_points(c1.y()-11.0) );
        write( diameter(points) );

        auto sf_tr = scaling_factor(points_tr);
        auto bar_tr = barycenter(points_tr);
        auto p0_tr = c2 + (points_tr[0] - bar_tr)/sf_tr;
        moveto( mm_to_points(p0_tr.x()), mm_to_points(p0_tr.y()) );
        for (size_t i = 1; i < points_tr.size(); i++) {
            auto v = (points_tr[i] - points_tr[i-1])/sf_tr;
            rlineto( mm_to_points(v.x()), mm_to_points(v.y()) );
        }
        closepath();

        moveto( mm_to_points(c2.x()-8.0), mm_to_points(c2.y()-11.0) );
        write( diameter(points_tr) );
        stroke();

        setrgbcolor(1.,0.,0.);
        linewidth(0.2);
        auto scale = std::max(m.col(0).norm(), m.col(1).norm());
        auto v0 = 8.0*m.col(0)/scale;
        auto v1 = 8.0*m.col(1)/scale;
        moveto( mm_to_points(c1.x()), mm_to_points(c1.y()) );
        rlineto( mm_to_points(v0(0)), mm_to_points(v0(1)) );
        closepath();
        moveto( mm_to_points(c1.x()), mm_to_points(c1.y()) );
        rlineto( mm_to_points(v1(0)), mm_to_points(v1(1)) );
        closepath();
        stroke();

        curr_index++;

        if (curr_index >= rows*cols) {
            curr_index = 0;
            show_page();
        }
    }

    ~postscript_plotter()
    {
        if (curr_index != 0) {
            stroke();
            show_page();
        }
    }
};

template<typename Mesh>
int
plot_inertia(Mesh& msh)
{
    using T = typename Mesh::coordinate_type;
    static const size_t DIM = Mesh::dimension;
    using point_type = typename Mesh::point_type;

    postscript_plotter psplot("inertia_transformations.ps");

    int elem_num = 0;
    for (auto& cl : msh)
    {
        Eigen::Matrix<T,DIM,DIM> m = scaled_inertia_axes(msh, cl);
        auto pts = points(msh, cl);
        std::vector<point_type> pp(pts.begin(), pts.end());
        psplot.add_poly(pp, m);
    }

    disk::silo_database db;
    db.create("inertia_starting_mesh.silo");
    db.add_mesh(msh, "mesh");

    return 0;
}

int main(int argc, char **argv)
{

    if (argc != 2)
    {
        std::cout << argv[0] << " <mesh file>" << std::endl;
        return 1;
    }

    auto mesh_filename = argv[1];

    return disk::dispatch_all_meshes(mesh_filename,
              [](auto ...args) { plot_inertia(args...); }
            );
}