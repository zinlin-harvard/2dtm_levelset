// CGAL must be installed for this to compile
// run cgal_create_CMakeLists -s gen_optimized_mesh
// cmake .      (might need -DCGAL_DIR=cgal_directory_path)
// make

#define CGAL_MESH_2_OPTIMIZER_VERBOSE
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/lloyd_optimize_mesh_2.h>
#include <iostream>
#include <string>
#include <sstream>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_mesh_vertex_base_2<K>                Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K>                  Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>        Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds>  CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT>            Criteria;
typedef CGAL::Delaunay_mesher_2<CDT, Criteria>              Mesher;
typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Point Point;

void add_cdt_constraints_from_file(const char *con_filename, CDT *cdt);
void write_cdt_to_gmsh(const CDT &cdt, const char *mesh_filename);
void add_x0(const std::list<double> &xy_pos, CDT *cdt, double x0);
void add_y0(const std::list<double> &xy_pos, CDT *cdt, double y0);
void add_x1(const std::list<double> &xy_pos, CDT *cdt, double x1);
void add_y1(const std::list<double> &xy_pos, CDT *cdt, double y1);
void generate_long_vec(std::list<double> &xy_pos, double xy0, double xy1, double el);
void generate_short_vec(std::list<double> &xy_pos, double xy0, double xy1);

int main(int argc, char **argv) {
    if (argc<4 || (argc>4 && argc<10)) {
        std::cerr << "Usage: gen_optimized_mesh  mesh_file_name  constraint_file_name  max_edge_length [bbx0 bby0 bbx1 bby1 per_x(0/1) per_y(0/1)] [min_angle (0.125)] [max_lloyd_iter (10)]" << std::endl;
        std::cerr << "    mesh output is .msh (gmsh) format" << std::endl;
        std::cerr << "    constraint file format is:" << std::endl;
        std::cerr << "        pointx pointy                   (for a point to include in mesh)" << std::endl;
        std::cerr << "        pointx pointy pointx2 pointy2   (for a line to included in mesh)" << std::endl;
        std::cerr << "    bb... refers to a \"bounding box\", which is a rectangle with lower" << std::endl;
        std::cerr << "        left corner bbx0,bby0 and upper right corner bbx1,bby1." << std::endl;
        std::cerr << "        per_x and per_y dictate whether the x and y boundaries have to be periodic" << std::endl;
        exit(EXIT_FAILURE);
        return -1;
    }

    const char *mesh_filename = argv[1];
    const char *con_filename = argv[2];
    double max_edge_length = atof(argv[3]);
    double min_angle = (argc>10) ? atof(argv[10]) : 0.125;
    unsigned max_lloyd_iter = (argc>11) ? atoi(argv[11]) : 10;
    std::cout << "max_lloyd_iter: " << max_lloyd_iter << std::endl;

    CDT cdt;
    add_cdt_constraints_from_file(con_filename, &cdt);
    std::cout << "Number of constrained vertices: " << cdt.number_of_vertices() << std::endl;
    std::cout << "Meshing..." << std::endl;
   
    if (argc>4) {
        double x0 = atof(argv[4]);
        double y0 = atof(argv[5]);
        double x1 = atof(argv[6]);
        double y1 = atof(argv[7]);
        if (x1<x0 || y1<y0) {
            std::cerr << "Must have x0<x1 and y0<y1. Exiting..." << std::endl;
            exit(EXIT_FAILURE);
        }
        bool per_x = (atoi(argv[8]) == 1) ? true : false;
        bool per_y = (atoi(argv[9]) == 1) ? true : false;

        std::list<double> x_vec, y_vec;
        if (per_x)
            generate_long_vec(x_vec, y0, y1, max_edge_length);
        else
            generate_short_vec(x_vec, y0, y1);
        if (per_y)
            generate_long_vec(y_vec, x0, x1, max_edge_length);
        else
            generate_short_vec(y_vec, x0, x1);
        
        add_x0(x_vec, &cdt, x0);
        add_y0(y_vec, &cdt, y0);
        add_x1(x_vec, &cdt, x1);
        add_y1(y_vec, &cdt, y1);
    }

    Mesher mesh(cdt);
    mesh.set_criteria(Criteria(min_angle, max_edge_length));
    mesh.refine_mesh();

    std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
    std::cout << "Run Lloyd optimization...";
    //CGAL::lloyd_optimize_mesh_2(cdt, CGAL::parameters::max_iteration_number = max_lloyd_iter);
    std::cout << " done." << std::endl;
    std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
    write_cdt_to_gmsh(cdt, mesh_filename);
}

void generate_long_vec(std::list<double> &xy_pos, double xy0, double xy1, double el) {
    double bnd_length = xy1 - xy0;
    double dxy_target = 0.5 * el;
    unsigned N = (int) floor(bnd_length / dxy_target); // N is # of line segments
    double dxy = bnd_length / double(N);
    std::cout << "dxy: " << dxy << std::endl;
    std::cout << "N: " << N << std::endl;

    for (unsigned cnt = 0; cnt<N+1; ++cnt)
        xy_pos.push_back(xy0 + cnt*dxy);
    xy_pos.back() = xy1; // last pos = xy1 (avoid additions that may not get xy1 exactly)
}

void generate_short_vec(std::list<double> &xy_pos, double xy0, double xy1) {
    xy_pos.push_back(xy0);
    xy_pos.push_back(xy1);
}

void add_x0(const std::list<double> &xy_pos, CDT *cdt, double x0) {
    std::list<double>::const_iterator xitc = xy_pos.begin();
    Vertex_handle v0, v1;
    while (xitc != --xy_pos.end()) { // x=0, increment y
        v0 = cdt->insert(Point(x0, *xitc));
        //std::cout << " 0 " << *xitc;
        v1 = cdt->insert(Point(x0, *(++xitc))); // update xitc
        //std::cout << " 0 " << *xitc << std::endl;
        cdt->insert_constraint(v0, v1);
    }
}

void add_y1(const std::list<double> &xy_pos, CDT *cdt, double y1) {
    std::list<double>::const_iterator xitc = xy_pos.begin();
    Vertex_handle v0, v1;
    while (xitc != --xy_pos.end()) { // x=0, increment y
        v0 = cdt->insert(Point(*xitc, y1));
        v1 = cdt->insert(Point(*(++xitc), y1)); // update xitc
        cdt->insert_constraint(v0, v1);
    }
}

void add_x1(const std::list<double> &xy_pos, CDT *cdt, double x1) {
    std::list<double>::const_reverse_iterator xitc = xy_pos.rbegin();
    Vertex_handle v0, v1;
    while (xitc != --xy_pos.rend()) { // x=0, increment y
        v0 = cdt->insert(Point(x1, *xitc));
        v1 = cdt->insert(Point(x1, *(++xitc))); // update xitc
        cdt->insert_constraint(v0, v1);
    }
}

void add_y0(const std::list<double> &xy_pos, CDT *cdt, double y0) {
    std::list<double>::const_reverse_iterator xitc = xy_pos.rbegin();
    Vertex_handle v0, v1;
    while (xitc != --xy_pos.rend()) { // x=0, increment y
        v0 = cdt->insert(Point(*xitc, y0));
        v1 = cdt->insert(Point(*(++xitc), y0)); // update xitc
        cdt->insert_constraint(v0, v1);
    }
}




// input file format:
//     point0x point0y                 (for a point to be added)
//     point0x point0y point1x point1y (for a line to be added)
//     ...
void add_cdt_constraints_from_file(const char *con_filename, CDT *cdt) {
    std::ifstream infile(con_filename);
    if (infile.is_open()) {
        double pointx0, pointy0, pointx1, pointy1;
        std::string line;
        while (std::getline(infile, line)) {
            std::stringstream linestream(line);
            linestream >> pointx0 >> pointy0;
            Vertex_handle v0 = cdt->insert(Point(pointx0, pointy0));
            if (linestream >> pointx1 >> pointy1) {
                Vertex_handle v1 = cdt->insert(Point(pointx1, pointy1));
                cdt->insert_constraint(v0, v1);
            } // else point already added by cdt->insert
        }
    }
}

void write_cdt_to_gmsh(const CDT &cdt, const char *mesh_filename) {
    // Follows the approach of output_surface_facets_to_off in 
    //   CGAL/IO/Complex_2_in_triangulation_3_file_writer.h
    std::ofstream outfile(mesh_filename);
    if (outfile.is_open()) {
        outfile << std::setprecision(16) 
            << "$MeshFormat \n2.2	0	8\n$EndMeshFormat \n";
        
        std::size_t nv = cdt.number_of_vertices();
        outfile << "$Nodes \n" << nv << " \n";
        std::map<Vertex_handle, int> Vi_map; // map each handle to a unique int
        unsigned i=1;                        // in gmsh, vertices start from 1
        CDT::Finite_vertices_iterator vit;
        for (vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit, ++i) {
            outfile << i << " " << vit->point() << " 0 \n"; // "0" for z=0
            Vi_map[vit] = i;
        }
        outfile << "$EndNodes\n";
        
        std::size_t nf = cdt.number_of_faces(); 
        outfile << "$Elements\n" << nf << "\n";
        CDT::Finite_faces_iterator fit;
        for (fit = cdt.finite_faces_begin(), i=1; fit != cdt.finite_faces_end(); ++fit, ++i)
            outfile << i << " 2 2 2 2"
                << " " << Vi_map[fit->vertex(0)]
                << " " << Vi_map[fit->vertex(1)]
                << " " << Vi_map[fit->vertex(2)] << "\n";
        outfile << "$EndElements";
    }
    // Outputs triangulation to "output.tri" (can check source in file_output methods in
    //   CGAL/Triangulation_data_structure_2.h). Was confused by numbering (which appears to
    //   have arisen from odd constraints originally imposed, so wrote above that outputs 
    //   directly to gmsh .msh format)
    //std::ofstream oFileT("output.tri", std::ios::out);
    //oFileT << std::setprecision(16) << cdt;
}

// old code to write out the points
/*    std::cout << "bnd_length: " << bnd_length << std::endl;*/
    //double x = 0, y = 0;
    //for (unsigned i=0; i<N; ++i) {
        //v0 = cdt.insert(Point(x, y));
        //v1 = cdt.insert(Point(x, y+dxy));
        //std::cout << x << " " << y << " " << x << " " << y+dxy << std::endl;
        //cdt.insert_constraint(v0, v1);
        //y += dxy;
    //}
    //for (unsigned i=0; i<N; ++i) {
        //v0 = cdt.insert(Point(x, y));
        //v1 = cdt.insert(Point(x+dxy, y));
        //std::cout << x << " " << y << " " << x+dxy << " " << y << std::endl;
        //cdt.insert_constraint(v0, v1);
        //x += dxy;
    //}
    //for (unsigned i=0; i<N; ++i) {
        //v0 = cdt.insert(Point(x, y));
        //v1 = cdt.insert(Point(x, y-dxy));
        //std::cout << x << " " << y << " " << x << " " << y-dxy << std::endl;
        //cdt.insert_constraint(v0, v1);
        //y -= dxy;
    //}
    //for (unsigned i=0; i<N-1; ++i) {
        //v0 = cdt.insert(Point(x, y));
        //v1 = cdt.insert(Point(x-dxy, y));
        //std::cout << x << " " << y << " " << x-dxy << " " << y << std::endl;
        //cdt.insert_constraint(v0, v1);
        //x -= dxy;
    /*}*/
    //v0 = cdt.insert(Point(x, y));
    //v1 = cdt.insert(Point(0, 0));
    //std::cout << x << " " << y << " " << 0 << " " << 0 << std::endl;
    //cdt.insert_constraint(v0, v1);
    /*double x = 0;*/
    //for (unsigned i=0; i<N; ++i) {
        //Vertex_handle v0 = cdt.insert(Point(x, y0));
        //Vertex_handle v1 = cdt.insert(Point(x+dxy, y0));
        //std::cout << x << " " << y0 << " " << x+dxy << " " << y0 << std::endl;
        //cdt.insert_constraint(v0, v1);
        //v0 = cdt.insert(Point(x, y1));
        //v1 = cdt.insert(Point(x+dxy, y1));
        //std::cout << x << " " << y1 << " " << x+dxy << " " << y1 << std::endl;
        //cdt.insert_constraint(v0, v1);
        //x += dxy;
    /*}*/


