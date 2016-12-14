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
#include <cstdio>

#include <unistd.h>
#include <limits.h>

/* By defining this the mesh is constructed in memory. This is needed when
 * this code has to be integrated somewhere. If the goal is to generate only
 * the meshfile, constructing in memory is useless, so data is written directly.
 */
//#define BUILD_IN_MEMORY

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
    printf("    -2: generate 2D mesh\n");
    printf("    -3: generate 3D mesh (default)\n");
    printf("    -x, -X: minimum and maximum x coordinate\n");
    printf("    -y, -Y: minimum and maximum y coordinate\n");
    printf("    -z, -Z: minimum and maximum z coordinate (ignored in 2D)\n");
    printf("    -m: number of elements in x direction\n");
    printf("    -n: number of elements in y direction\n");
    printf("    -p: number of elements in z direction (ignored in 2D)\n");
}

template<typename T, size_t DIM>
struct mesh_params;

template<typename T>
struct mesh_params<T, 2>
{
    mesh_params()
    {
        nx = 1; ny = 1;

        min_x = 0.0; max_x = 1.0;
        min_y = 0.0; max_y = 1.0;
    }

    char *      filename;
    size_t      nx, ny;
    T           min_x, max_x, min_y, max_y;
};

template<typename T>
struct mesh_params<T, 3>
{
    mesh_params()
    {
        nx = 1; ny = 1; nz = 1;
        min_x = 0.0; max_x = 1.0;
        min_y = 0.0; max_y = 1.0;
        min_z = 0.0; max_z = 1.0;
    }

    char *      filename;
    size_t      nx, ny, nz;
    T           min_x, max_x, min_y, max_y, min_z, max_z;
};

template<typename T, size_t DIM>
struct mesh_data;

template<typename T>
struct mesh_data<T, 2>
{
    typedef size_t                      idx_type;
    typedef std::array<T, 2>            point_type;
    typedef std::array<idx_type, 4>     cell_type;
    typedef std::array<idx_type, 2>     bnd_type;

#ifdef BUILD_IN_MEMORY
    std::vector<point_type>             points;
    std::vector<cell_type>              cells;
    std::vector<bnd_type>               boundaries;
#endif
};

template<typename T>
struct mesh_data<T, 3>
{
    typedef size_t                      idx_type;
    typedef std::array<T, 3>            point_type;
    typedef std::array<idx_type, 8>     cell_type;
    typedef std::array<idx_type, 4>     bnd_type;

#ifdef BUILD_IN_MEMORY
    std::vector<point_type>             points;
    std::vector<cell_type>              cells;
    std::vector<bnd_type>               boundaries;
#endif
};

template<typename T>
#ifdef BUILD_IN_MEMORY
mesh_data<T,2>
#else
bool
#endif
generate(const mesh_params<T,2>& mp)
{
#ifdef BUILD_IN_MEMORY
    mesh_data<T,2> md;
    return md;
#else
    return true;
#endif
}

#ifdef BUILD_IN_MEMORY
template<typename T, size_t DIM>
void
write_points(const mesh_data<T,DIM>& md, FILE *fp)
{
    fprintf(fp, "%lu\n", md.points.size());
    for (auto& p : md.points)
    {
        if (DIM == 2) fprintf(fp, "%.15f %.15f\n", p[0], p[1]);
        if (DIM == 3) fprintf(fp, "%.15f %.15f %.15f\n", p[0], p[1], p[2]);
    }
}

template<typename T, size_t DIM>
void
write_cells(const mesh_data<T,DIM>& md, FILE *fp)
{
    fprintf(fp, "%lu\n", md.cells.size());
    for (auto& c : md.cells)
    {
        if (DIM == 2) fprintf(fp, "%lu %lu %lu %lu\n", c[0], c[1], c[2], c[3]);
        if (DIM == 3) fprintf(fp, "%lu %lu %lu %lu %lu %lu %lu %lu\n",
                                  c[0], c[1], c[2], c[3],
                                  c[4], c[5], c[6], c[7]);
    }
}

template<typename T, size_t DIM>
void
write_boundaries(const mesh_data<T,DIM>& md, FILE *fp)
{
    fprintf(fp, "%lu\n", md.boundaries.size());
    for (auto& b : md.boundaries)
    {
        if (DIM == 2) fprintf(fp, "%lu %lu\n", b[0], b[1]);
        if (DIM == 3) fprintf(fp, "%lu %lu %lu %lu\n", b[0], b[1], b[2], b[3]);
    }
}

template<typename T, size_t DIM>
bool
write_mesh(const mesh_params<T, DIM>& mp, const mesh_data<T, DIM>& md)
{
    FILE *fp = stdout;

    if (mp.filename) fp = fopen(mp.filename, "w");

    if (!fp)
    {
        std::cout << "Unable to open " << mp.filename << std::endl;
        return false;
    }

    write_points(md, fp);
    write_cells(md, fp);
    write_boundaries(md, fp);

    if (mp.filename) fclose(fp);
    return true;
}
#endif

template<typename T>
#ifdef BUILD_IN_MEMORY
mesh_data<T,3>
#else
bool
#endif
generate(const mesh_params<T,3>& mp)
{
#ifdef BUILD_IN_MEMORY
    mesh_data<T,3> md;
#else
    FILE *fp = stdout;
    if (mp.filename) fp = fopen(mp.filename, "w");
    if (!fp)
    {
        std::cout << "Unable to open " << mp.filename << std::endl;
        return false;
    }
#endif

    auto hx = (mp.max_x - mp.min_x)/mp.nx;
    auto hy = (mp.max_y - mp.min_y)/mp.ny;
    auto hz = (mp.max_z - mp.min_z)/mp.nz;

    auto num_points = (mp.nx+1)*(mp.ny+1)*(mp.nz+1);
    auto num_vols = mp.nx * mp.ny * mp.nz;
    auto num_bf = 2*mp.nx*mp.ny + 2*mp.nx*mp.nz + 2*mp.ny*mp.nz;

#ifdef BUILD_IN_MEMORY
    md.points.reserve(num_points);
    md.cells.reserve(num_vols);
    md.boundaries.reserve(num_bf);
#else
    fprintf(fp, "%lu\n", num_points);
#endif

    /* Generate points */
    for (size_t k = 0; k < mp.nz+1; k++)
    {
        for (size_t j = 0; j < mp.ny+1; j++)
        {
            for (size_t i = 0; i < mp.nx+1; i++)
            {
                typename mesh_data<T,3>::point_type p;
                p[0] = mp.min_x + i*hx;
                p[1] = mp.min_y + j*hy;
                p[2] = mp.min_z + k*hz;

#ifdef BUILD_IN_MEMORY
                md.points.push_back(p);
#else
                fprintf(fp, "%.15f %.15f %.15f\n", p[0], p[1], p[2]);
#endif
            }
        }
    }

    /* Generate volumes */
    size_t xy_offset = (mp.nx+1)*(mp.ny+1);
    size_t y_offset = mp.ny+1;

#ifndef BUILD_IN_MEMORY
    fprintf(fp, "%lu\n", num_vols);
#endif
    for (size_t k = 0; k < mp.nz; k++)
    {
        for (size_t j = 0; j < mp.ny; j++)
        {
            for (size_t i = 0; i < mp.nx; i++)
            {
                typename mesh_data<T,3>::cell_type   vol;
                vol[0] =   k   * xy_offset +   j   * y_offset +   i;
                vol[1] =   k   * xy_offset +   j   * y_offset + (i+1);
                vol[2] =   k   * xy_offset + (j+1) * y_offset +   i;
                vol[3] =   k   * xy_offset + (j+1) * y_offset + (i+1);
                vol[4] = (k+1) * xy_offset +   j   * y_offset +   i;
                vol[5] = (k+1) * xy_offset +   j   * y_offset + (i+1);
                vol[6] = (k+1) * xy_offset + (j+1) * y_offset +   i;
                vol[7] = (k+1) * xy_offset + (j+1) * y_offset + (i+1);

#ifdef BUILD_IN_MEMORY
                md.cells.push_back(vol);
#else
                fprintf(fp, "%lu %lu %lu %lu %lu %lu %lu %lu\n",
                            vol[0], vol[1], vol[2], vol[3],
                            vol[4], vol[5], vol[6], vol[7]);
#endif
            }
        }
    }

#ifndef BUILD_IN_MEMORY
    fprintf(fp, "%lu\n", num_bf);
#endif
    /* Generate boundary faces */
    for (size_t k = 0; k < mp.nz; k++)
    {
        for (size_t j = 0; j < mp.ny; j++)
        {
            for (size_t i = 0; i < mp.nx; i++)
            {
                if (i == 0)
                {
                    typename mesh_data<T,3>::bnd_type   bf;
                    bf[0] =   k   * xy_offset +   j   * y_offset +   i;
                    bf[1] =   k   * xy_offset + (j+1) * y_offset +   i;
                    bf[2] = (k+1) * xy_offset + (j+1) * y_offset +   i;
                    bf[3] = (k+1) * xy_offset +   j   * y_offset +   i;
#ifdef BUILD_IN_MEMORY
                    md.boundaries.push_back(bf);
#else
#endif
                }

                if (i == mp.nx-1)
                {
                    typename mesh_data<T,3>::bnd_type   bf;
                    bf[0] =   k   * xy_offset +   j   * y_offset + (i+1);
                    bf[1] =   k   * xy_offset + (j+1) * y_offset + (i+1);
                    bf[2] = (k+1) * xy_offset + (j+1) * y_offset + (i+1);
                    bf[3] = (k+1) * xy_offset +   j   * y_offset + (i+1);
#ifdef BUILD_IN_MEMORY
                    md.boundaries.push_back(bf);
#else
                    fprintf(fp, "%lu %lu %lu %lu\n", bf[0], bf[1], bf[2], bf[3]);
#endif
                }

                if (j == 0)
                {
                    typename mesh_data<T,3>::bnd_type   bf;
                    bf[0] =   k   * xy_offset +   j   * y_offset +   i;
                    bf[1] = (k+1) * xy_offset +   j   * y_offset +   i;
                    bf[2] = (k+1) * xy_offset +   j   * y_offset + (i+1);
                    bf[3] =   k   * xy_offset +   j   * y_offset + (i+1);
#ifdef BUILD_IN_MEMORY
                    md.boundaries.push_back(bf);
#else
                    fprintf(fp, "%lu %lu %lu %lu\n", bf[0], bf[1], bf[2], bf[3]);
#endif
                }

                if (j == mp.ny-1)
                {
                    typename mesh_data<T,3>::bnd_type   bf;
                    bf[0] =   k   * xy_offset + (j+1) * y_offset +   i;
                    bf[1] = (k+1) * xy_offset + (j+1) * y_offset +   i;
                    bf[2] = (k+1) * xy_offset + (j+1) * y_offset + (i+1);
                    bf[3] =   k   * xy_offset + (j+1) * y_offset + (i+1);
#ifdef BUILD_IN_MEMORY
                    md.boundaries.push_back(bf);
#else
                    fprintf(fp, "%lu %lu %lu %lu\n", bf[0], bf[1], bf[2], bf[3]);
#endif
                }

                if (k == 0)
                {
                    typename mesh_data<T,3>::bnd_type   bf;
                    bf[0] =   k   * xy_offset +   j   * y_offset +   i;
                    bf[1] =   k   * xy_offset +   j   * y_offset + (i+1);
                    bf[2] =   k   * xy_offset + (j+1) * y_offset + (i+1);
                    bf[3] =   k   * xy_offset + (j+1) * y_offset +   i;
#ifdef BUILD_IN_MEMORY
                    md.boundaries.push_back(bf);
#else
                    fprintf(fp, "%lu %lu %lu %lu\n", bf[0], bf[1], bf[2], bf[3]);
#endif
                }

                if (k == mp.nz-1)
                {
                    typename mesh_data<T,3>::bnd_type   bf;
                    bf[0] = (k+1) * xy_offset +   j   * y_offset +   i;
                    bf[1] = (k+1) * xy_offset +   j   * y_offset + (i+1);
                    bf[2] = (k+1) * xy_offset + (j+1) * y_offset + (i+1);
                    bf[3] = (k+1) * xy_offset + (j+1) * y_offset +   i;
#ifdef BUILD_IN_MEMORY
                    md.boundaries.push_back(bf);
#else
                    fprintf(fp, "%lu %lu %lu %lu\n", bf[0], bf[1], bf[2], bf[3]);
#endif
                }
            }
        }
    }
#ifdef BUILD_IN_MEMORY
    return md;
#else
    if (mp.filename) fclose(fp);
    return true;
#endif
}

int main(int argc, char **argv)
{
    using T = double;

    size_t      nx = 1, ny = 1, nz = 1;
    T           min_x = 0.0, max_x = 1.0;
    T           min_y = 0.0, max_y = 1.0;
    T           min_z = 0.0, max_z = 1.0;
    bool        three_dimensions = true;
    int         ch;

    while ( (ch = getopt(argc, argv, "x:X:y:Y:z:Z:m:n:p:23")) != -1 )
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
            case '2':   three_dimensions = false;
            case '3':   three_dimensions = true;

            case 'h':
            case '?':
            default:
                usage(argv[0]);
                exit(1);
        }
    }

    argc -= optind;
    argv += optind;

    if (three_dimensions)
    {
        mesh_params<T,3> mp;
        mp.min_x = min_x; mp.max_x = max_x;
        mp.min_y = min_y; mp.max_x = max_y;
        mp.min_z = min_z; mp.max_x = max_z;
        mp.nx = nx; mp.ny = ny; mp.nz = nz;

        if ( mp.max_x < mp.min_x or mp.max_y < mp.min_y or mp.max_z < mp.min_z )
        {
            std::cout << "Invalid range" << std::endl;
            return 1;
        }

        mp.filename = argv[0];

#ifdef BUILD_IN_MEMORY
        mesh_data<T,3> md = generate(mp);
        write_mesh(mp, md);
#else
        generate(mp);
#endif
    }
    else
    {
        mesh_params<T,2> mp;
        mp.min_x = min_x; mp.max_x = max_x;
        mp.min_y = min_y; mp.max_x = max_y;
        mp.nx = nx; mp.ny = ny;

        if ( mp.max_x < mp.min_x or mp.max_y < mp.min_y )
        {
            std::cout << "Invalid range" << std::endl;
            return 1;
        }

        mp.filename = argv[0];

#ifdef BUILD_IN_MEMORY
        mesh_data<T,2> md = generate(mp);
        write_mesh(mp, md);
#else
        generate(mp);
#endif
    }

    return 0;
}
