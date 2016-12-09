/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016 - matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code for scientific publications, you are required to
 * cite it.
 */


/*
 * This is a trivial cartesian mesh generator. Nothing fancy, but does the job.
 * The output is Netgen-like and INDICES START FROM ZERO:
 *
 *      <number of points>
 *      x1 y1 z1
 *      x2 y2 z2
 *      ...
 *      xn yn zn
 *      <number of cells>
 *      node1 node2 node3 ... node8
 *      <number of boundary faces>
 *      node1 node2 ... node4
 *
 * In the hexahedron the nodes are named using a 3-bit binary number: 0 if
 * the coordinate is not present, 1 if it is (convert in binary numbers in
 * figure and look at the name of the axis).
 *
 *
 *          ^ Y
 *          |                           Face 0: 0 2 6 4
 *        2 _____________ 3             Face 1: 1 3 7 5
 *         /            /|              Face 2: 0 1 3 2
 *      6 / |        7 / |              Face 3: 4 5 7 6
 *       /____________/  |              Face 4: 0 4 5 1
 *      |   |         |  |      X       Face 5: 2 6 7 3
 *      | 0 _ _ _ _ _ |_ |1 ---->
 *      |  /          |  /
 *      |             | /
 *      |/____________|/
 *     4              5
 *     /
 *    v Z
 *
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <array>
#include <cstdlib>

#include <unistd.h>
#include <limits.h>

template<typename T>
void
output_mesh(const std::vector<std::array<T, 3>>& points,
            const std::vector<std::array<size_t, 8>>& volumes,
            const std::vector<std::array<size_t, 4>>& boundary_faces,
            std::ostream& os)
{
    os << points.size() << std::endl;
    for (auto& p : points)
    {
        os << std::setprecision(15) << std::fixed;
        for (auto& c : p)
            os << c << "\t";
        os << std::endl;
    }

    os << volumes.size() << std::endl;
    for (auto& v : volumes)
    {
        for (auto& c : v)
            os << c << " ";
        os << std::endl;
    }

    os << boundary_faces.size() << std::endl;
    for (auto& bf : boundary_faces)
    {
        for (auto& c : bf)
            os << c << " ";
        os << std::endl;
    }
}

template<typename T>
T strto(const char * nptr, char ** endptr = nullptr);

template<>
float strto(const char * nptr, char ** endptr)
{
    return strtof(nptr, endptr);
}

template<>
double strto(const char * nptr, char ** endptr)
{
    return strtod(nptr, endptr);
}

template<>
unsigned long strto(const char * nptr, char ** endptr)
{
    return strtoul(nptr, endptr, 10);
}

void
usage(const char *progname)
{
    printf("Usage: %s <options> <filename>\n\n", progname);
    printf("    -x, -X: minimum and maximum x coordinate\n");
    printf("    -y, -Y: minimum and maximum y coordinate\n");
    printf("    -z, -Z: minimum and maximum z coordinate\n");
    printf("    -m: number of elements in x direction\n");
    printf("    -n: number of elements in y direction\n");
    printf("    -p: number of elements in z direction\n");
}

int main(int argc, char **argv)
{
    using T = double;

    size_t nx = 1;
    size_t ny = 1;
    size_t nz = 1;

    T   min_x = 0.0, max_x = 1.0,
        min_y = 0.0, max_y = 1.0,
        min_z = 0.0, max_z = 1.0;

    int ch;

    while ( (ch = getopt(argc, argv, "x:X:y:Y:z:Z:m:n:p:")) != -1 )
    {
        switch(ch)
        {
            case 'm':   nx = strto<size_t>(optarg); break;
            case 'n':   ny = strto<size_t>(optarg); break;
            case 'p':   nz = strto<size_t>(optarg); break;
            case 'x':   min_x = strto<T>(optarg); break;
            case 'X':   max_x = strto<T>(optarg); break;
            case 'y':   min_y = strto<T>(optarg); break;
            case 'Y':   max_y = strto<T>(optarg); break;
            case 'z':   min_z = strto<T>(optarg); break;
            case 'Z':   max_z = strto<T>(optarg); break;

            case 'h':
            case '?':
            default:
                usage(argv[0]);
                exit(1);
        }
    }

    argc -= optind;
    argv += optind;

    if ( max_x < min_x or max_y < min_y or max_z < min_z )
    {
        std::cout << "Invalid range" << std::endl;
        return 1;
    }

    char *filename = argv[0];

    std::vector<std::array<T, 3>>           points;
    std::vector<std::array<size_t, 8>>      volumes;
    std::vector<std::array<size_t, 4>>      boundary_faces;

    double hx = (max_x - min_x)/nx;
    double hy = (max_y - min_y)/ny;
    double hz = (max_z - min_z)/nz;

    auto num_points = (nx+1)*(ny+1)*(nz+1);
    auto num_vols = nx*ny*nz;
    auto num_bf = 2*nx*ny + 2*nx*nz + 2*ny*nz;

    points.reserve(num_points);
    volumes.reserve(num_vols);
    boundary_faces.reserve(num_bf);

    /* Generate points */
    for (size_t k = 0; k < nz+1; k++)
    {
        for (size_t j = 0; j < ny+1; j++)
        {
            for (size_t i = 0; i < nx+1; i++)
            {
                std::array<T, 3> p;
                p[0] = min_x + i*hx;
                p[1] = min_y + j*hy;
                p[2] = min_z + k*hz;

                points.push_back(p);
            }
        }
    }

    /* Generate volumes */
    size_t xy_offset = (nx+1)*(ny+1);
    size_t y_offset = ny+1;

    for (size_t k = 0; k < nz; k++)
    {
        for (size_t j = 0; j < ny; j++)
        {
            for (size_t i = 0; i < nx; i++)
            {
                std::array<size_t, 8>   vol;
                vol[0] =   k   * xy_offset +   j   * y_offset +   i;
                vol[1] =   k   * xy_offset +   j   * y_offset + (i+1);
                vol[2] =   k   * xy_offset + (j+1) * y_offset +   i;
                vol[3] =   k   * xy_offset + (j+1) * y_offset + (i+1);
                vol[4] = (k+1) * xy_offset +   j   * y_offset +   i;
                vol[5] = (k+1) * xy_offset +   j   * y_offset + (i+1);
                vol[6] = (k+1) * xy_offset + (j+1) * y_offset +   i;
                vol[7] = (k+1) * xy_offset + (j+1) * y_offset + (i+1);

                volumes.push_back(vol);
            }
        }
    }

    /* Generate boundary faces */
    for (size_t k = 0; k < nz; k++)
    {
        for (size_t j = 0; j < ny; j++)
        {
            for (size_t i = 0; i < nx; i++)
            {
                if (i == 0)
                {
                    std::array<size_t, 4>   bf;
                    bf[0] =   k   * xy_offset +   j   * y_offset +   i;
                    bf[1] =   k   * xy_offset + (j+1) * y_offset +   i;
                    bf[2] = (k+1) * xy_offset + (j+1) * y_offset +   i;
                    bf[3] = (k+1) * xy_offset +   j   * y_offset +   i;
                    boundary_faces.push_back(bf);
                }

                if (i == nx-1)
                {
                    std::array<size_t, 4>   bf;
                    bf[0] =   k   * xy_offset +   j   * y_offset + (i+1);
                    bf[1] =   k   * xy_offset + (j+1) * y_offset + (i+1);
                    bf[2] = (k+1) * xy_offset + (j+1) * y_offset + (i+1);
                    bf[3] = (k+1) * xy_offset +   j   * y_offset + (i+1);
                    boundary_faces.push_back(bf);
                }

                if (j == 0)
                {
                    std::array<size_t, 4>   bf;
                    bf[0] =   k   * xy_offset +   j   * y_offset +   i;
                    bf[1] = (k+1) * xy_offset +   j   * y_offset +   i;
                    bf[2] = (k+1) * xy_offset +   j   * y_offset + (i+1);
                    bf[3] =   k   * xy_offset +   j   * y_offset + (i+1);
                    boundary_faces.push_back(bf);
                }

                if (j == ny-1)
                {
                    std::array<size_t, 4>   bf;
                    bf[0] =   k   * xy_offset + (j+1) * y_offset +   i;
                    bf[1] = (k+1) * xy_offset + (j+1) * y_offset +   i;
                    bf[2] = (k+1) * xy_offset + (j+1) * y_offset + (i+1);
                    bf[3] =   k   * xy_offset + (j+1) * y_offset + (i+1);
                    boundary_faces.push_back(bf);
                }

                if (k == 0)
                {
                    std::array<size_t, 4>   bf;
                    bf[0] =   k   * xy_offset +   j   * y_offset +   i;
                    bf[1] =   k   * xy_offset +   j   * y_offset + (i+1);
                    bf[2] =   k   * xy_offset + (j+1) * y_offset + (i+1);
                    bf[3] =   k   * xy_offset + (j+1) * y_offset +   i;
                    boundary_faces.push_back(bf);
                }

                if (k == nz-1)
                {
                    std::array<size_t, 4>   bf;
                    bf[0] = (k+1) * xy_offset +   j   * y_offset +   i;
                    bf[1] = (k+1) * xy_offset +   j   * y_offset + (i+1);
                    bf[2] = (k+1) * xy_offset + (j+1) * y_offset + (i+1);
                    bf[3] = (k+1) * xy_offset + (j+1) * y_offset +   i;
                    boundary_faces.push_back(bf);
                }
            }
        }
    }

    if (filename)
    {
        std::ofstream ofs(filename);
        if (!ofs.is_open())
        {
            std::cout << "Unable to open " << filename << std::endl;
            exit(1);
        }

        output_mesh(points, volumes, boundary_faces, ofs);
        ofs.close();
    }
    else
    {
        output_mesh(points, volumes, boundary_faces, std::cout);
    }

    return 0;
}
