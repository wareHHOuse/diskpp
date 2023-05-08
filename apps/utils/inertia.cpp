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

    template<typename T, size_t DIM>
    T diameter(const std::vector<disk::point<T,DIM>>& points)
    {
        T diam = 0.0;
        for (size_t i = 0; i < points.size(); i++)
            for (size_t j = i; j < points.size(); j++)
                diam = std::max( diam, distance(points[i], points[j]) );
        
        return diam;
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

    void make_grid(void)
    {
        ps_ofs << "0.2 setlinewidth" << std::endl;

        for (size_t row = 0; row < rows+1; row++)
        {
            double line_start_x = mm_to_points(margin_x);
            double line_start_y = mm_to_points(margin_y + row*box_height);

            double line_end_x = mm_to_points(margin_x + cols*box_width);
            double line_end_y = mm_to_points(margin_y + row*box_height);

            ps_ofs << line_start_x << " " << line_start_y << " moveto" << std::endl;
            ps_ofs << line_end_x << " " << line_end_y << " lineto" << std::endl;
        }

        for (size_t col = 0; col < cols+1; col++)
        {
            double line_start_x = mm_to_points(margin_x + col*box_width);
            double line_start_y = mm_to_points(margin_y);

            double line_end_x = mm_to_points(margin_x + col*box_width);
            double line_end_y = mm_to_points(margin_y + rows*box_height);

            ps_ofs << line_start_x << " " << line_start_y << " moveto" << std::endl;
            ps_ofs << line_end_x << " " << line_end_y << " lineto" << std::endl;
        }
    }

    void stroke(void)
    {
        ps_ofs << "stroke" << std::endl;
    }

    void show_page(void)
    {
        ps_ofs << "showpage" << std::endl;
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
    void add_poly(const std::vector<disk::point<T,3>>& points, const std::vector<disk::point<T,3>>&)
    {}

    template<typename T>
    void add_poly(const std::vector<disk::point<T,2>>& points,
        const std::vector<disk::point<T,2>>& points_tr)
    {
        if (curr_index == 0) {
            make_grid();
            stroke();
            ps_ofs << "1 setlinewidth" << std::endl;
        }
        
        auto [c1, c2] = centers<T>(curr_index);
        
        auto sf_orig = scaling_factor(points);
        auto bar_orig = barycenter(points);
        auto p0 = c1 + (points[0] - bar_orig)/sf_orig;
        ps_ofs << mm_to_points(p0.x()) << " " << mm_to_points(p0.y()) << " moveto" << std::endl;
        for (size_t i = 1; i < points.size(); i++) {
            auto v = (points[i] - points[i-1])/sf_orig;
            ps_ofs << mm_to_points(v.x()) << " " << mm_to_points(v.y()) << " rlineto" << std::endl;
        }
        ps_ofs << "closepath" << std::endl;
        ps_ofs << mm_to_points(c1.x()-10.0) << " " << mm_to_points(c1.y()-10.0) << " moveto" << std::endl;
        ps_ofs << "(" << diameter(points) << ") show" << std::endl;

        auto sf_tr = scaling_factor(points_tr);
        auto bar_tr = barycenter(points_tr);
        auto p0_tr = c2 + (points_tr[0] - bar_tr)/sf_tr;
        ps_ofs << mm_to_points(p0_tr.x()) << " " << mm_to_points(p0_tr.y()) << " moveto" << std::endl;
        for (size_t i = 1; i < points_tr.size(); i++) {
            auto v = (points_tr[i] - points_tr[i-1])/sf_tr;
            ps_ofs << mm_to_points(v.x()) << " " << mm_to_points(v.y()) << " rlineto" << std::endl;
        }
        ps_ofs << "closepath" << std::endl;
        ps_ofs << "stroke" << std::endl;
        ps_ofs << mm_to_points(c2.x()-10.0) << " " << mm_to_points(c2.y()-10.0) << " moveto" << std::endl;
        ps_ofs << "(" << diameter(points_tr) << ") show" << std::endl;

        curr_index++;

        if (curr_index >= rows*cols) {
            curr_index = 0;
            stroke();
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




template<typename T>
void
plot_ps(const std::vector<T>& points, const T& bar, const std::string& filename)
{
    std::ofstream ofs(filename);

    double pt_scale = 2.54/72.;

    ofs << "newpath" << std::endl;
    ofs << "144 144 moveto" << std::endl;

    auto p0 = points[0] - bar;

    ofs << int(10.*p0.x()/pt_scale) << " " << int(10.*p0.y()/pt_scale) << " rmoveto" << std::endl;

    for (size_t i = 1; i < points.size(); i++) {
        const auto& pt = (points[i] - points[i-1]);
        int x = 10.*pt.x()/pt_scale;
        int y = 10.*pt.y()/pt_scale;
        ofs << x << " " << y << " rlineto" << std::endl;
    }
    ofs << "closepath" << std::endl;
    //ofs << "8 setlinewidth" << std::endl;
    ofs << "stroke" << std::endl;
    ofs << "showpage" << std::endl;
}


template<typename Mesh>
int
plot_inertia(Mesh& msh)
{
    using T = typename Mesh::coordinate_type;
    static const size_t DIM = Mesh::dimension;
    using point_type = typename Mesh::point_type;

    std::vector<double> axis1_x, axis1_y, axis2_x, axis2_y;
    std::vector<double> elem_nums;

    postscript_plotter psplot("test.ps");

    int elem_num = 0;
    for (auto& cl : msh)
    {
        Eigen::Matrix<T,DIM,DIM> m = scaled_inertia_axes(msh, cl);
        axis1_x.push_back( m(0,0) );
        axis1_y.push_back( m(1,0) );
        axis2_x.push_back( m(0,1) );
        axis2_y.push_back( m(1,1) );
        elem_nums.push_back( elem_num++ );
        auto pts = points(msh, cl);

        auto bar = barycenter(msh, cl);

        std::vector<point_type> pp(pts.begin(), pts.end());
        //plot_ps(pp, bar, "elem_" + std::to_string(elem_num) + ".ps" );

        auto tr = [&](const point_type& pt) {
            Eigen::Matrix<T, DIM, 1> trp = m.transpose()*pt.to_vector();
            point_type ret;
            ret.x() = trp(0);
            ret.y() = trp(1);
            return ret;
        };

        std::vector<point_type> pt(pp.size());
        std::transform(pp.begin(), pp.end(), pt.begin(), tr);
        psplot.add_poly(pp, pt);
    }

    disk::silo_database db;
    db.create("inertia.silo");

    db.add_mesh(msh, "mesh");

    disk::silo_zonal_variable<double> a1x("a1x", axis1_x);
    disk::silo_zonal_variable<double> a1y("a1y", axis1_y);
    disk::silo_zonal_variable<double> a2x("a2x", axis2_x);
    disk::silo_zonal_variable<double> a2y("a2y", axis2_y);
    disk::silo_zonal_variable<double> elnum("elnum", elem_nums);

    db.add_variable("mesh", a1x);
    db.add_variable("mesh", a1y);
    db.add_variable("mesh", a2x);
    db.add_variable("mesh", a2y);
    db.add_variable("mesh", elnum);

    db.add_expression("a1", "{a1x, a1y}", DB_VARTYPE_VECTOR);
    db.add_expression("a2", "{a2x, a2y}", DB_VARTYPE_VECTOR);

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