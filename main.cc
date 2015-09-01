#include <exception>
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <set>
#include <vector>

#include <ctime>
#include <cmath>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
//#include <boost/thread.hpp>

//#include <CGAL/convex_hull_2.h>
#include <CGAL/exceptions.h>
//#include <CGAL/Arrangement_on_surface_2.h>

#include "visibility.h"
#include "errors.h"
#include "model.h"

using std::list;
using std::set;
using std::string;
using std::vector;

//const int threads_num = 4;
//boost::mutex coutmutex;
//boost::mutex totalsummutex;

/*
void do_sum(const list<ProbPolygon>& prob_polys1,
            const list<ProbPolygon>& prob_polys2,
            const Obstacles& obs, mpfr_class *total_visb) {

//    for (list<ProbPolygon>::const_iterator it = prob_polys1.begin(); it != prob_polys1.end(); ++it)
//        for (list<ProbPolygon>::const_iterator jt = prob_polys2.begin(); jt != prob_polys2.end(); ++jt) {
//            coutmutex.lock();
//            std::cout << std::setprecision(10) << std::endl
//                      << "Polygon 1: " << it->poly << std::endl
//                      << "Polygon 2: " << jt->poly << std::endl;
//            coutmutex.unlock();
//        }

    mpfr_class partial_sum(0);
    for (list<ProbPolygon>::const_iterator it = prob_polys1.begin(); it != prob_polys1.end(); ++it)
        for (list<ProbPolygon>::const_iterator jt = prob_polys2.begin(); jt != prob_polys2.end(); ++jt) {
            try {
                coutmutex.lock();
                std::cout << std::setprecision(10) << std::endl
                          << "Polygon 1: " << it->poly << std::endl
                          << "Polygon 2: " << jt->poly << std::endl;
                coutmutex.unlock();
                mpfr_class visb = volume(it->poly, jt->poly, obs);
                partial_sum += visb*it->weight*jt->weight;
            } catch (exception_base& e) {
                coutmutex.lock();
                std::cerr << "Error!" << std::endl;
                if( error_num const * c = boost::get_error_info<errno_code>(e) )
                    std::cerr << "OS says: " << error_to_str(*c) << std::endl;

                if( std::string const * ds = boost::get_error_info<err_description>(e) )
                    std::cerr << *ds ;

                if( int const * cl = boost::get_error_info<code_line>(e) )
                    std::cerr << " line " << *cl << std::endl;

                std::cerr << std::setprecision(10) << std::endl
                << "Polygon 1: " << it->poly << std::endl
                << "Polygon 2: " << jt->poly << std::endl;
                coutmutex.unlock();
                exit(EXIT_FAILURE);
            } catch (CGAL::Failure_exception& e) {
                coutmutex.lock();
                std::cerr << "CGAL Error: " << e.what() << std::endl;
                std::cerr << std::setprecision(10) << std::endl
                << "Polygon 1: " << it->poly << std::endl
                << "Polygon 2: " << jt->poly << std::endl;
                coutmutex.unlock();
                exit(EXIT_FAILURE);
            } catch (std::exception& e) {
                coutmutex.lock();
                std::cerr << "Exception: " << e.what() << std::endl;
                std::cerr << std::setprecision(10) << std::endl
                << "Polygon 1: " << it->poly << std::endl
                << "Polygon 2: " << jt->poly << std::endl;
                coutmutex.unlock();
                exit(EXIT_FAILURE);
            }
        }

    *total_visb += partial_sum;
}

int run_monkeys_multithread(int argc, char *argv[]) {
    double error;
    if ( argc != 2 ) {
        std::cout<< "usage: " << argv[0] << " error\n";
        return 0;
    } else {
        error = std::atof(argv[1]);
    }

    time_t rawtime;
    time (&rawtime);

    string file_name;
    boost::filesystem::path dir(std::to_string(rawtime));
    if (boost::filesystem::create_directory(dir))
        file_name = std::to_string(rawtime) + "/result";
    else
        file_name = string("result") + std::to_string(rawtime);

    std::ofstream f1(file_name + ".csv");
    f1 << "#id, sigma, epsilon, visibility, time" << std::endl;
    f1.close();

    double sigma = 10;
    Model model(error);
    try {
        model.read_obstacles_from_file("../data/obstacles1.txt");
        model.set_first_monkey(-40, 0, sigma);
        model.set_second_monkey(40, 0, sigma);
    } catch (exception_base& e) {
        std::cerr << "Error!" << std::endl;
        if( error_num const * c = boost::get_error_info<errno_code>(e) )
            std::cerr << "OS says: " << error_to_str(*c) << std::endl;

        if( std::string const * ds = boost::get_error_info<err_description>(e) )
            std::cerr << *ds ;

        if( int const * cl = boost::get_error_info<code_line>(e) )
            std::cerr << " line " << *cl << std::endl;
    } catch (CGAL::Failure_exception& e) {
        std::cerr << "CGAL Error: " << e.what() << std::endl;
    } catch (std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
    }

    list<ProbPolygon> monkey1 = model.monkey(0);
    list<ProbPolygon> monkey2 = model.monkey(1);

    list<SimplePolygon> obs_polys;
    for (list<Obstacle>::const_iterator it = model.obstacles().obstacles().begin();
         it != model.obstacles().obstacles().end(); ++it)
        obs_polys.push_back(it->obstacle);

    print_to_file(-40, 0, 40, 0, sigma, monkey1, monkey2, obs_polys, file_name + ".svg");
    std::ofstream of1(file_name + "_polygons1.csv");
    for (list<ProbPolygon>::iterator it = monkey1.begin(); it != monkey1.end(); ++it)
        of1 << it->poly << ", " << it->weight << std::endl;
    std::ofstream of2(file_name + "_polygons2.csv");
    of1.close();
    for (list<ProbPolygon>::iterator it = monkey2.begin(); it != monkey2.end(); ++it)
        of2 << it->poly << ", " << it->weight << std::endl;
    of2.close();

    model.trim_obstacles();

    vector<mpfr_class*> part_sums;
    boost::thread threads[threads_num];
    vector<list<ProbPolygon>* > monkey1_parts;
    vector<list<ProbPolygon>* > monkey2_parts;
    vector<Obstacles> obstacles_parts;

    for (int i = 0; i < threads_num; ++i)
        part_sums.push_back(new mpfr_class(0));

    for (int i = 0; i < threads_num; ++i) {
        monkey1_parts.push_back(new list<ProbPolygon>());
        monkey2_parts.push_back(new list<ProbPolygon>());
        monkey2_parts.back()->insert(monkey2_parts.back()->begin(), monkey2.begin(), monkey2.end());
        obstacles_parts.push_back(Obstacles());
    }

    for (struct {list<ProbPolygon>::iterator it; int i;} loop = {monkey1.begin(), 0}; loop.it != monkey1.end(); ++(loop.it), ++(loop.i)) {
        ProbPolygon prob_poly;
        prob_poly.weight = loop.it->weight;
        for (SimplePolygon::Vertex_const_iterator v_it = loop.it->poly.vertices_begin(); v_it != loop.it->poly.vertices_end(); ++v_it)
            prob_poly.poly.push_back(Point(CGAL::to_double(v_it->x()), CGAL::to_double(v_it->y())));
        monkey1_parts[loop.i%threads_num]->push_back(prob_poly);
    }

    for (list<ProbPolygon>::iterator it = monkey2.begin(); it != monkey2.end(); ++it) {
        ProbPolygon prob_poly;
        prob_poly.weight = it->weight;
        for (SimplePolygon::Vertex_const_iterator v_it = it->poly.vertices_begin(); v_it != it->poly.vertices_end(); ++v_it)
            prob_poly.poly.push_back(Point(CGAL::to_double(v_it->x()), CGAL::to_double(v_it->y())));
        for (int i = 0; i < threads_num; ++i)
            monkey2_parts[i]->push_back(prob_poly);
    }

    for (list<Obstacle>::const_iterator it = model.obstacles().obstacles().begin(); it != model.obstacles().obstacles().end(); ++it) {
        SimplePolygon poly;
        for (SimplePolygon::Vertex_const_iterator v_it = it->obstacle.vertices_begin(); v_it != it->obstacle.vertices_end(); ++v_it)
            poly.push_back(Point(CGAL::to_double(v_it->x()), CGAL::to_double(v_it->y())));
        for (int i = 0; i < threads_num; ++i)
            obstacles_parts[i].addObstacle(poly);
    }

    for (int i = 0; i < threads_num; ++i)
        threads[i]= boost::thread(do_sum, *(monkey1_parts[i]), *(monkey2_parts[i]), obstacles_parts[i], part_sums[i]);

    for (int i = 0; i < threads_num; ++i)
        threads[i].join();

    mpfr_class total_visb(0);

    for (int i = 0; i < threads_num; ++i) {
        total_visb += *part_sums[i];
    }

    for (int i = 0; i < threads_num; ++i) {
        delete part_sums[i];
        delete monkey1_parts[i];
        delete monkey2_parts[i];
    }

    time_t rawtime_end;
    time (&rawtime_end);

    f1.open(file_name + ".csv", std::ofstream::out | std::ofstream::app | std::ofstream::ate);
    f1 << rawtime << ", " << sigma << ", " << error << ", " << total_visb.get_d() << ", " << rawtime_end-rawtime << std::endl;
    f1.close();
    
    return 0;
}
*/
/*
void generate_obstacles() {
    list<SimplePolygon> obs;
    for (int j = 0; j < 300; ++j) {
        SimplePolygon ob;
        double x = 0.1*rand___.get_int(-500, 500);
        double y = 0.1*rand___.get_int(-500, 500);
        double r = rand___.get_double(0.5, 0.75);
        double angle = rand___.get_double(0,M_PI_2);
        for (int i = 0; i < 5; ++i) {
            ob.push_back(Point(x+r*sin(angle+0.4*M_PI*i),y+r*cos(angle+0.4*M_PI*i)));
        }
        obs.push_back(ob);
    }

    Model model(0);
    model.set_obstacles(obs);
    obs.clear();


    std::ofstream f2("../data/obstacles.txt");
    for (list<Obstacle>::const_iterator it = model.obstacles().obstacles().begin();
         it != model.obstacles().obstacles().end(); ++it) {
        f2 << it->obstacle;
        obs.push_back(it->obstacle);
    }
    f2.close();

    print_to_file(SimplePolygon(), SimplePolygon(), obs, "../data/obstacles.svg");
}
*/
/*
void run_monkeys_complete(int argc, char *argv[]) {
    double error;
    string obstacles_filename;
    double x1, y1, x2, y2;
    double sigma;
    if ( argc != 8 ) {
        std::cout<< "usage: " << argv[0] << " error obstacles_filename x1 y1 x2 y2 sigma\n";
        return;
    } else {
        error = std::atof(argv[1]);
        obstacles_filename = string(argv[2]);
        x1 = std::atof(argv[3]);
        y1 = std::atof(argv[4]);
        x2 = std::atof(argv[5]);
        y2 = std::atof(argv[6]);
        sigma = std::atof(argv[7]);
    }

    time_t rawtime;
    time (&rawtime);

    string file_name;
    boost::filesystem::path dir(std::to_string(rawtime));
    if (boost::filesystem::create_directory(dir))
        file_name = std::to_string(rawtime) + "/result";
    else
        file_name = string("result") + std::to_string(rawtime);

    std::ofstream f1(file_name + ".csv");
    f1 << "#id, sigma, epsilon, visibility, time" << std::endl;
    f1.close();

    mpfr_class total_visb(0);
    Model model(error);
    try {
//        model.set_obstacles(obs);
        model.read_obstacles_from_file(obstacles_filename);
        model.set_first_monkey(x1, y1, sigma);
        model.set_second_monkey(x2, y2, sigma);
    } catch (exception_base& e) {
        std::cerr << "Error!" << std::endl;
        if( error_num const * c = boost::get_error_info<errno_code>(e) )
            std::cerr << "OS says: " << error_to_str(*c) << std::endl;

        if( std::string const * ds = boost::get_error_info<err_description>(e) )
            std::cerr << *ds ;

        if( int const * cl = boost::get_error_info<code_line>(e) )
            std::cerr << " line " << *cl << std::endl;
    } catch (CGAL::Failure_exception& e) {
        std::cerr << "CGAL Error: " << e.what() << std::endl;
    } catch (std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
    }

    list<ProbPolygon> monkey1 = model.monkey(0);
    list<ProbPolygon> monkey2 = model.monkey(1);

    list<SimplePolygon> obs_polys;
    for (list<Obstacle>::const_iterator it = model.obstacles().obstacles().begin();
         it != model.obstacles().obstacles().end(); ++it)
        obs_polys.push_back(it->obstacle);

    double opacity = 1-pow(1-0.2, 2*(exp(error)-1));
    print_to_file(x1, y1, x2, y2, sigma, monkey1, monkey2, opacity, obs_polys, file_name + ".svg");
    std::ofstream of_log(file_name + "_log.txt");
    of_log << "monkey1=(" << x1 << ", " << y1 << ")" << std::endl
        << "monkey2=(" << x2 << ", " << y2 << ")" << std::endl
        << "sigma=" << sigma << std::endl
        << "epsilon=" << error << std::endl;
    of_log.close();
    std::ofstream of1(file_name + "_polygons1.csv");
    for (list<ProbPolygon>::iterator it = monkey1.begin(); it != monkey1.end(); ++it)
        of1 << it->poly << ", " << it->weight << std::endl;
    std::ofstream of2(file_name + "_polygons2.csv");
    of1.close();
    for (list<ProbPolygon>::iterator it = monkey2.begin(); it != monkey2.end(); ++it)
        of2 << it->poly << ", " << it->weight << std::endl;
    of2.close();

    model.trim_obstacles();

//    unsigned long total_count = monkey1.size()*monkey2.size();
//    int count = 0;
    for (list<ProbPolygon>::iterator it = monkey1.begin(); it != monkey1.end(); ++it)
        for (list<ProbPolygon>::iterator jt = monkey2.begin(); jt != monkey2.end(); ++jt) {
            ProbPolygon& prob_poly1 = *it;
            ProbPolygon& prob_poly2 = *jt;
            try {
                mpfr_class visb = volume(prob_poly1.poly, prob_poly2.poly, model.obstacles());
                total_visb += visb*prob_poly1.weight*prob_poly2.weight;
//                std::cout << "\r" << ++count << " of " << total_count << std::flush;
//                double vol4 = CGAL::to_double(prob_poly1.poly.area()*prob_poly2.poly.area());
//                std::cout << "Visibility = " << visb/vol4 << std::endl << std::endl;
            } catch (exception_base& e) {
                std::ofstream f_err(file_name + "_exception.txt");

                std::cerr << "Error!" << std::endl;
                if( error_num const * c = boost::get_error_info<errno_code>(e) ) {
                    std::cerr << "OS says: " << error_to_str(*c) << std::endl;
                    f_err << error_to_str(*c) << std::endl;
                }

                if( std::string const * ds = boost::get_error_info<err_description>(e) ) {
                    std::cerr << *ds ;
                    f_err << *ds;
                }

                if( int const * cl = boost::get_error_info<code_line>(e) ) {
                    std::cerr << " line " << *cl << std::endl;
                    f_err << " line " << *cl << std::endl;
                }

                std::cerr << "Polygons written to " << file_name << "_exception_polygons.txt"
                            << std::setprecision(10) << std::endl
                            << "Polygon 1: " << prob_poly1.poly << std::endl
                            << "Polygon 2: " << prob_poly2.poly << std::endl;

                f_err << prob_poly1.poly << prob_poly2.poly;
                f_err.close();

                return;
            } catch (CGAL::Failure_exception& e) {
                std::cerr << "CGAL Error: " << e.what() << std::endl;
                std::cerr << "Polygons written to " << file_name << "_exception_polygons.txt"
                            << std::setprecision(10) << std::endl
                            << "Polygon 1: " << prob_poly1.poly << std::endl
                            << "Polygon 2: " << prob_poly2.poly << std::endl;

                std::ofstream f_err(file_name + "_exception.txt");
                f_err << e.what() << std::endl;
                f_err << prob_poly1.poly << prob_poly2.poly;
                f_err.close();

                return;
            } catch (std::exception& e) {
                std::cerr << "Exception: " << e.what() << std::endl;
                std::cerr << "Polygons written to " << file_name << "_exception.txt"
                          << std::setprecision(10) << std::endl
                          << "Polygon 1: " << prob_poly1.poly << std::endl
                          << "Polygon 2: " << prob_poly2.poly << std::endl;

                std::ofstream f_err(file_name + "_exception.txt");
                f_err << e.what() << std::endl;
                f_err << prob_poly1.poly << prob_poly2.poly;
                f_err.close();

                return;
            }
        }
//    std::cout << std::endl;

    time_t rawtime_end;
    time (&rawtime_end);

    f1.open(file_name + ".csv", std::ofstream::out | std::ofstream::app | std::ofstream::ate);
    f1 << rawtime << ", " << sigma << ", " << error << ", " << total_visb.get_d() << ", " << rawtime_end-rawtime << std::endl;
    f1.close();

//        std::cout << std::endl << "Visibility = " << total_visb << std::endl << std::endl;

}
*/
/*
void generate_monkeys(int argc, char *argv[]) {
    double error;
    if ( argc != 2 ) {
        std::cout<< "usage: " << argv[0] << " error\n";
        return;
    } else {
        error = std::atof(argv[1]);
    }

//    time_t rawtime;
//    time (&rawtime);

    string file_name1[4];
    string file_name2[4];
    for (int i = 0; i < 4; ++i) {
        file_name1[i] = "monkey1_part" + std::to_string(i) + ".txt";
        file_name2[i] = "monkey2_part" + std::to_string(i) + ".txt";
    }
//    boost::filesystem::path dir(std::to_string(rawtime));
//    if (boost::filesystem::create_directory(dir))
//        file_name = std::to_string(rawtime) + "/result";
//    else
//        file_name = string("result") + std::to_string(rawtime);

//    std::ofstream f1(file_name + ".csv");
//    f1 << "#id, sigma, epsilon, visibility, time" << std::endl;
//    f1.close();

//    mpfr_class total_visb(0);
    double sigma = 10;
    Model model(error);
    try {
        model.read_obstacles_from_file("../data/obstacles1.txt");
        //        model.set_obstacles(obs);
        model.set_first_monkey(-40, 0, sigma);
        model.set_second_monkey(40, 0, sigma);
    } catch (exception_base& e) {
        std::cerr << "Error!" << std::endl;
        if( error_num const * c = boost::get_error_info<errno_code>(e) )
            std::cerr << "OS says: " << error_to_str(*c) << std::endl;

        if( std::string const * ds = boost::get_error_info<err_description>(e) )
            std::cerr << *ds ;

        if( int const * cl = boost::get_error_info<code_line>(e) )
            std::cerr << " line " << *cl << std::endl;
    } catch (CGAL::Failure_exception& e) {
        std::cerr << "CGAL Error: " << e.what() << std::endl;
    } catch (std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
    }

    list<ProbPolygon> monkey1 = model.monkey(0);
    list<ProbPolygon> monkey2 = model.monkey(1);

    list<SimplePolygon> obs_polys;
    for (list<Obstacle>::const_iterator it = model.obstacles().obstacles().begin();
         it != model.obstacles().obstacles().end(); ++it)
        obs_polys.push_back(it->obstacle);

    print_to_file(-40, 0, 40, 0, sigma, monkey1, monkey2, obs_polys, "scenario.svg");

    std::ofstream of1[4];
    std::ofstream of2[4];
    for (int i = 0; i < 4; ++i) {
        of1[i].open(file_name1[i]);
        of2[i].open(file_name2[i]);
    }

    int i = 0;
    for (list<ProbPolygon>::iterator it = monkey1.begin(); it != monkey1.end(); ++it)
        of1[i++%4] << it->poly << it->weight << std::endl;

    i = 0;
    for (list<ProbPolygon>::iterator it = monkey2.begin(); it != monkey2.end(); ++it)
        of2[i++%4] << it->poly << it->weight << std::endl;

    for (int i = 0; i < 4; ++i) {
        of1[i].close();
        of2[i].close();
    }
}

void run_monkeys(int argc, char *argv[]) {
    int i1, i2;
    if ( argc != 3 ) {
        std::cout<< "usage: " << argv[0] << " num1 num2\n";
        return;
    } else {
        i1 = std::atoi(argv[1]);
        i2 = std::atoi(argv[2]);
    }

    time_t rawtime_begin;
    time (&rawtime_begin);

    std::ifstream f_obs("../data/obstacles1.txt");
    SimplePolygon ob;

    Obstacles obs;
    while (f_obs >> ob) {
        obs.addObstacle(ob);
    }
    f_obs.close();

    list<ProbPolygon> monkey1;
    list<ProbPolygon> monkey2;

    ProbPolygon prob_poly;

    std::ifstream f1("monkey1_part" + std::to_string(i1) + ".txt");
    while (f1 >> prob_poly.poly >> prob_poly.weight) {
        monkey1.push_back(prob_poly);
        prob_poly.poly.clear();
    }
    f1.close();

    std::ifstream f2("monkey2_part" + std::to_string(i2) + ".txt");
    while (f2 >> prob_poly.poly >> prob_poly.weight) {
        monkey2.push_back(prob_poly);
        prob_poly.poly.clear();
    }
    f2.close();

    //    unsigned long total_count = monkey1.size()*monkey2.size();
    //    int count = 0;
    mpfr_class total_visb(0);
    std::ofstream f_out("result_"+std::to_string(i1)+"_"+std::to_string(i2)+".csv");
    f_out << "#i1, i2, sigma,  Visibility" << std::endl;
    f_out.close();
    for (list<ProbPolygon>::iterator it = monkey1.begin(); it != monkey1.end(); ++it)
        for (list<ProbPolygon>::iterator jt = monkey2.begin(); jt != monkey2.end(); ++jt) {
            ProbPolygon& prob_poly1 = *it;
            ProbPolygon& prob_poly2 = *jt;
//            std::cout << prob_poly1.poly << " \t\t\t " << prob_poly1.weight << std::endl;
//            std::cout << prob_poly2.poly << " \t\t\t " << prob_poly2.weight << std::endl;
            try {
                mpfr_class visb = volume(prob_poly1.poly, prob_poly2.poly, obs);
                total_visb += visb*prob_poly1.weight*prob_poly2.weight;
                //                std::cout << "\r" << ++count << " of " << total_count << std::flush;
                //                double vol4 = CGAL::to_double(prob_poly1.poly.area()*prob_poly2.poly.area());
                //                std::cout << "Visibility = " << visb/vol4 << std::endl << std::endl;
            } catch (exception_base& e) {
                std::cerr << "Error!" << std::endl;
                if( error_num const * c = boost::get_error_info<errno_code>(e) )
                    std::cerr << "OS says: " << error_to_str(*c) << std::endl;

                if( std::string const * ds = boost::get_error_info<err_description>(e) )
                    std::cerr << *ds ;

                if( int const * cl = boost::get_error_info<code_line>(e) )
                    std::cerr << " line " << *cl << std::endl;

                std::cerr << std::setprecision(10) << std::endl
                << "Polygon 1: " << prob_poly1.poly << std::endl
                << "Polygon 2: " << prob_poly2.poly << std::endl;

                std::ofstream f_err("exception_polygons.txt");
                //                CGAL::set_binary_mode(f_err);
                f_err << prob_poly1.poly << prob_poly2.poly;
                f_err.close();

                return;
            } catch (CGAL::Failure_exception& e) {
                std::cerr << "CGAL Error: " << e.what() << std::endl;
                std::cerr << std::setprecision(10) << std::endl
                << "Polygon 1: " << prob_poly1.poly << std::endl
                << "Polygon 2: " << prob_poly2.poly << std::endl;

                std::ofstream f_err("exception_polygons.txt");
                //                CGAL::set_binary_mode(f_err);
                f_err << prob_poly1.poly << prob_poly2.poly;
                f_err.close();

                return;
            } catch (std::exception& e) {
                std::cerr << "Exception: " << e.what() << std::endl;
                std::cerr << std::setprecision(10) << std::endl
                << "Polygon 1: " << prob_poly1.poly << std::endl
                << "Polygon 2: " << prob_poly2.poly << std::endl;

                std::ofstream f_err("exception_polygons.txt");
                //                CGAL::set_binary_mode(f_err);
                f_err << prob_poly1.poly << prob_poly2.poly;
                f_err.close();

                return;
            }
        }
    //    std::cout << std::endl;

    time_t rawtime_end;
    time (&rawtime_end);

    f_out.open("result_"+std::to_string(i1)+"_"+std::to_string(i2)+".csv", std::ofstream::out | std::ofstream::app | std::ofstream::ate);
    f_out << i1 << ", " << i2 << ", " << total_visb.get_d() << ", " << rawtime_end-rawtime_begin << std::endl;
    f_out.close();
    
    //        std::cout << std::endl << "Visibility = " << total_visb << std::endl << std::endl;
}
*/
/*
int run_random_polygons(int argc, char *argv[]) {
    for (int i = 0; i < 100; ++i) {
        list<Point> pts1, pts2;
        set<double> xs;

        while (pts1.size() < 5 || pts2.size() < 5) {
            double x1 = rand___.get_int(-100, 100);
            double y1 = rand___.get_int(-100, 100);

            if (xs.find(x1) == xs.end()) {
                xs.insert(x1);
                pts1.push_back(Point(x1, y1));
            }

            double x2 = rand___.get_int(-100, 100);
            double y2 = rand___.get_int(-100, 100);
            if (xs.find(x2) == xs.end()) {
                xs.insert(x2);
                pts2.push_back(Point(x2, y2));
            }
        }

        SimplePolygon polygon1, polygon2;
        CGAL::convex_hull_2(pts1.begin(), pts1.end(), std::back_inserter(polygon1));
        CGAL::convex_hull_2(pts2.begin(), pts2.end(), std::back_inserter(polygon2));

        std::ofstream f1("data/data.txt");
        f1 << polygon1 << polygon2;
        f1.close();

        Obstacles obs;
        list<SimplePolygon> obs_polys;
        std::ofstream f2("data/obstacles.txt");
        for (int j = 0; j < 30; ++j) {
            SimplePolygon ob;
            double x = 0.5*rand___.get_int(-100, 100);
            if (xs.find(x) != xs.end() || xs.find(x-1) != xs.end() || xs.find(x+1) != xs.end())
                continue;

            xs.insert(x);
            xs.insert(x+1);
            xs.insert(x-1);

            double y = rand___.get_int(-50, 50);
            ob.push_back(Point(x+0,y+0));
            ob.push_back(Point(x+1,y+2));
            ob.push_back(Point(x-1,y+2));
            obs.addObstacle(ob);
            obs_polys.push_back(ob);
            f2 << ob;
        }
        f2.close();

        list<SimplePolygon> monkey[2];

        PolygonSet poly_set[2];
        poly_set[0].insert(polygon1);
        poly_set[1].insert(polygon2);
        for (int i = 0; i < 2; ++i) {
            for (list<Obstacle>::const_iterator it = obs.obstacles().begin(); it != obs.obstacles().end(); ++it) {
                poly_set[i].difference(it->obstacle);
            }
            list<Polygon> polys;
            poly_set[i].polygons_with_holes(std::back_inserter(polys));
            for (list<Polygon>::iterator it = polys.begin(); it != polys.end(); ++it) {
                list<SimplePolygon> convex_polys;
                partition(*it, &convex_polys);
                monkey[i].insert(monkey[i].end(), convex_polys.begin(), convex_polys.end());
            }
        }

        print_to_file(monkey[0], monkey[1], obs_polys, "aaa" + std::to_string(i) +".svg");

        mpfr_class total_visb(0);
        for (SimplePolygon poly1 : monkey[0])
            for (SimplePolygon poly2 : monkey[1]) {
                try {
                    total_visb += volume(poly1, poly2, obs);
                } catch (exception_base& e) {
                    std::cerr << "Error!" << std::endl;
                    if( error_num const * c = boost::get_error_info<errno_code>(e) )
                        std::cerr << "OS says: " << error_to_str(*c) << std::endl;

                    if( std::string const * ds = boost::get_error_info<err_description>(e) )
                        std::cerr << *ds ;

                    if( int const * cl = boost::get_error_info<code_line>(e) )
                        std::cerr << " line " << *cl << std::endl;

                    std::cerr << std::setprecision(10) << std::endl
                    << "Polygon 1: " << poly1 << std::endl
                    << "Polygon 2: " << poly2 << std::endl;

                } catch (CGAL::Failure_exception& e) {
                    std::cerr << "CGAL Error: " << e.what() << std::endl;
                    std::cerr << std::setprecision(10) << std::endl
                    << "Polygon 1: " << poly1 << std::endl
                    << "Polygon 2: " << poly2 << std::endl;
                } catch (std::exception& e) {
                    std::cerr << "Exception: " << e.what() << std::endl;
                    std::cerr << std::setprecision(10) << std::endl
                    << "Polygon 1: " << poly1 << std::endl
                    << "Polygon 2: " << poly2 << std::endl;
                }
            }

        std::cout << "Visibility = " << (total_visb/mpfr_class(polygon1.area()*polygon2.area())).get_d() << std::endl;
    }

    return 0;
}

void run_polygons(int argc, char *argv[]) {
    SimplePolygon polygon1, polygon2;
    std::ifstream f1("../data/data.txt");
    f1 >> polygon1 >> polygon2;
    f1.close();

    std::cout << "Polygon 1: " << polygon1 << std::endl
              << "Polygon 2: " << polygon2 << std::endl;

    Model model(0);
    model.read_obstacles_from_file("../data/obstacles1.txt");
    model.write_obstacles_to_file("../data/obstacles_temp.txt");

    std::cout << model.obstacles().obstacles().size() << std::endl;

    list<SimplePolygon> monkey[2];

    PolygonSet poly_set[2];
    poly_set[0].insert(polygon1);
    poly_set[1].insert(polygon2);
    for (int i = 0; i < 2; ++i) {
        for (list<Obstacle>::const_iterator it = model.obstacles().obstacles().begin(); it != model.obstacles().obstacles().end(); ++it) {
            poly_set[i].difference(it->obstacle);
        }
        list<Polygon> polys;
        poly_set[i].polygons_with_holes(std::back_inserter(polys));
        for (list<Polygon>::iterator it = polys.begin(); it != polys.end(); ++it) {
            list<SimplePolygon> convex_polys;
            partition(*it, &convex_polys);
            monkey[i].insert(monkey[i].end(),convex_polys.begin(),convex_polys.end());
        }
    }

    list<SimplePolygon> obs_polys;
    for (list<Obstacle>::const_iterator it = model.obstacles().obstacles().begin(); it != model.obstacles().obstacles().end(); ++it) {
        obs_polys.push_back(it->obstacle);
    }
    print_to_file(monkey[0], monkey[1], obs_polys, "aaa.svg");

    mpfr_class total_visb(0);
    for (SimplePolygon poly1 : monkey[0])
        for (SimplePolygon poly2 : monkey[1])
            try {
                total_visb += volume(poly1, poly2, model.obstacles());
            } catch (exception_base& e) {
                std::cerr << "Error!" << std::endl;
                if( error_num const * c = boost::get_error_info<errno_code>(e) )
                    std::cerr << "OS says: " << error_to_str(*c) << std::endl;

                if( std::string const * ds = boost::get_error_info<err_description>(e) )
                    std::cerr << *ds ;

                if( int const * cl = boost::get_error_info<code_line>(e) )
                    std::cerr << " line " << *cl << std::endl;

                std::cerr << std::setprecision(10) << std::endl
                << "Polygon 1: " << poly1 << std::endl
                << "Polygon 2: " << poly2 << std::endl;

                return;
            } catch (CGAL::Failure_exception& e) {
                std::cerr << "CGAL Error: " << e.what() << std::endl;
                std::cerr << std::setprecision(10) << std::endl
                << "Polygon 1: " << poly1 << std::endl
                << "Polygon 2: " << poly2 << std::endl;

                return;
            } catch (std::exception& e) {
                std::cerr << "Exception: " << e.what() << std::endl;
                std::cerr << std::setprecision(10) << std::endl
                << "Polygon 1: " << poly1 << std::endl
                << "Polygon 2: " << poly2 << std::endl;

                return;
            }

    std::cout << "Visibility = " << total_visb.get_d() << std::endl;
    std::cout << "Area = " << polygon1.area()*polygon2.area() << std::endl;

}


void run_segments(int argc, char *argv[]) {
    Obstacles obs;
    list<SimplePolygon> obs_polys;
    std::ifstream f2("../data/obstacles_temp.txt");
    SimplePolygon ob;

    while (f2 >> ob) {
        obs.addObstacle(ob);
        obs_polys.push_back(ob);
    }

    std::cout << obs.obstacles().size() << std::endl;

//    print_to_file(monkey[0], monkey[1], obs_polys, "aaa.svg");

    try {
        Segment segs[4];// = { Segment(Point(-29.5, -2), Point(-24.47, -8.98)),
//            Segment(Point(-24.47, -8.98), Point(-28.5, 0)),
//            Segment(Point(30.63, 12.04), Point(30.03981281, 11)),
//            Segment(Point(33.13613042, 12.72773916), Point(30.63, 12.04)) };
        std::stringstream ss;
        ss << "-29.5 -2 -24.47 -8.98 -24.47 -8.98 -28.5 0 30.63 12.04 30.03981281 11 33.13613042 12.72773916 30.63 12.04";

//        int foo,bar;
        ss >> segs[0] >> segs[1] >> segs[2] >> segs[3];
//        segs[0] << "-29.5 -2 -24.47 -8.98";
//        segs[1] << "-24.47 -8.98 -28.5 0";
//        segs[2] << "30.63 12.04 30.03981281 11";
//        segs[3] << "33.13613042 12.72773916 30.63 12.04";

        SimplePolygon p1, p2;
        p1.push_back(segs[0].source());
        p1.push_back(segs[0].target());
        p1.push_back(segs[1].target());
        p1.push_back(segs[1].source());

        p2.push_back(segs[2].source());
        p2.push_back(segs[2].target());
        p2.push_back(segs[3].target());
        p2.push_back(segs[3].source());

        print_to_file(p1, p2, obs_polys, "aaa.svg");

        for (int i = 0; i < 4; ++i)
            std::cout << "[" << segs[i].source().x() << "," << segs[i].source().y() << "];" << std::endl
                      << "[" << segs[i].target().x() << "," << segs[i].target().y() << "];" << std::endl;

        std::cout << (volume(segs, obs)).get_d() << std::endl;
    } catch (exception_base& e) {
        std::cerr << "Error!" << std::endl;
        if( error_num const * c = boost::get_error_info<errno_code>(e) )
            std::cerr << "OS says: " << error_to_str(*c) << std::endl;

        if( std::string const * ds = boost::get_error_info<err_description>(e) )
            std::cerr << *ds ;

        if( int const * cl = boost::get_error_info<code_line>(e) )
            std::cerr << " line " << *cl << std::endl;

        return;
    } catch (CGAL::Failure_exception& e) {
        std::cerr << "CGAL Error: " << e.what() << std::endl;

        return;
    } catch (std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;

        return;
    }
}
*/
/*
int main() {
    InexactSegment s1(InexactPoint(0, -31), InexactPoint(11.3636, -38.5116));
    InexactSegment s2(InexactPoint(11.3636, -38.5116), InexactPoint(12.2357, -37.6247));
    InexactSegment s3(InexactPoint(20.642, -43.4444), InexactPoint(12.2357, -37.6247));
    InexactSegment s4(InexactPoint(20.642, -43.4444), InexactPoint(63, 14));
    InexactSegment segs[4] = { s1, s2, s3, s4 };
    InexactLine l1(InexactPoint(0, -31), InexactPoint(12.2357, -37.6247));
    InexactLine l2(InexactPoint(11.3636, -38.5116), InexactPoint(20.642, -43.4444));
    InexactLine l3(InexactPoint(0, -31), InexactPoint(20.642, -43.4444));
    InexactLine l4(InexactPoint(11.3636, -38.5116), InexactPoint(63, 14));

    std::cout << visibility_intl::intersects_in_order(l1, segs) << std::endl;
    std::cout << visibility_intl::intersects_in_order(l2, segs) << std::endl;
    std::cout << visibility_intl::intersects_in_order(l3, segs) << std::endl;
    std::cout << visibility_intl::intersects_in_order(l4, segs) << std::endl;
}
*/

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

void run_point_sample(int argc, char *argv[]) {
    string obstacles_filename;
    double x1, y1, x2, y2;
    //double sigma = 15;
    double sigma = 3;
    if ( argc != 6 ) {
        std::cout<< "usage: " << argv[0] << "obstacles_filename x1 y1 x2 y2\n";
        return;
    } else {
        obstacles_filename = string(argv[1]);
        x1 = std::atof(argv[2]);
        y1 = std::atof(argv[3]);
        x2 = std::atof(argv[4]);
        y2 = std::atof(argv[5]);
    }

    time_t rawtime;
    time (&rawtime);

    string file_name;
    boost::filesystem::path dir(boost::to_string(rawtime));
    if (boost::filesystem::create_directory(dir))
        file_name = boost::to_string(rawtime) + "/result";
    else
        file_name = string("result") + boost::to_string(rawtime);

    std::ofstream of_log(file_name + "_log.txt");
    of_log << "monkey1=(" << x1 << ", " << y1 << ")" << std::endl
            << "monkey2=(" << x2 << ", " << y2 << ")" << std::endl
            << "sigma=" << sigma << std::endl;
    of_log.close();

    Model model(0);
    model.read_obstacles_from_file(obstacles_filename);

    typedef boost::normal_distribution<double> NormalDistribution;
    typedef boost::mt19937 RandomGenerator;
    typedef boost::variate_generator<RandomGenerator&, NormalDistribution> GaussianGenerator;

    static RandomGenerator rng(static_cast<unsigned> (time(0)));
    NormalDistribution norm1(0, sigma);
    NormalDistribution norm2(0, sigma);

    GaussianGenerator generator1(rng, norm1);
    GaussianGenerator generator2(rng, norm2);

    std::ofstream flog(file_name+"_point_sample_method.txt");
    flog << "num pairs, prob" << std::endl;
    int see_sum = 0;
    for (int i = 0; i <= 100000; ++i) {
        double rnd_x1 = x1 + generator1();
        double rnd_y1 = y1 + generator1();
        double rnd_x2 = x2 + generator2();
        double rnd_y2 = y2 + generator2();

        Segment seg(Point(rnd_x1, rnd_y1), Point(rnd_x2, rnd_y2));

        bool see = true;
        for (list<Obstacle>::const_iterator it = model.obstacles().obstacles().begin(); it != model.obstacles().obstacles().end(); ++it) {
            if (!(it->obstacle.bounded_side(seg.source()) == CGAL::ON_UNBOUNDED_SIDE &&
                it->obstacle.bounded_side(seg.target()) == CGAL::ON_UNBOUNDED_SIDE)) {
                see = false;
                break;
            }

            for (SimplePolygon::Edge_const_iterator jt = it->obstacle.edges_begin(); jt != it->obstacle.edges_end(); ++jt) {
                if (CGAL::do_intersect(seg, *jt)) {
                    see = false;
                    break;
                }
            }

            if (see == false)
                break;
        }

        if (see) {
            see_sum++;
        }

        flog << i << ", " << double(see_sum)/double(i+1) << std::endl;
    }

    flog.close();
}

int main(int argc, char *argv[]) {
//    generate_obstacles();
//    run_polygons(argc, argv);
//    run_segments(argc, argv);
//    run_monkeys_complete(argc, argv);
//    run_point_sample(argc,argv);
//    run_monkeys_multithread(argc, argv);
//    generate_monkeys(argc, argv);
//    run_monkeys(argc, argv);

    std::map<int,double> vals;
    for (double error = 0.0001; error <= 0.5; error += 0.0001 ) {
        Model model(error);
        model.set_first_monkey(0, 0, 15);
        list<ProbPolygon> polys = model.monkey(0);
        int size = 0;
        for (ProbPolygon p : polys) {
            size += p.poly.size();
        }

        if (vals.find(size) == vals.end()) {
            vals[size] = error;
            std::cout << size << "," << error << std::endl;
        }
//        std::cout << error << " " << size << std::endl;
    }


//    std::ifstream f("data/utrecht_parsed.txt");
//    SimplePolygon ob;
////    Bbox b(225, 450, 400, 470);
//    int size = 0;
//    while (f >> ob) {
////        b += ob.bbox();
//      size += ob.size();
//    }
//
//    std::cout << size << std::endl;
/*
    std::ifstream f("../data/utrecht.txt");
    string line;
    list<SimplePolygon> obs;
    list<Point> pts;
    while (! f.eof() ) {
        //place the line from myfile into the
        //line variable:
        std::getline(f, line);

        if (line.empty()) {
            SimplePolygon ob(pts.begin(), pts.end());

            if (ob.is_clockwise_oriented())
                ob.reverse_orientation();

            list<CGAL::Partition_traits_2<Kernel>::Polygon_2> partition_polys;
            CGAL::optimal_convex_partition_2(ob.vertices_begin(),
                                             ob.vertices_end(),
                                             std::back_inserter(partition_polys));
            if (!CGAL::partition_is_valid_2(ob.vertices_begin(),
                                            ob.vertices_end(),
                                            partition_polys.begin(),
                                            partition_polys.end()))
                throw boost::enable_current_exception(geometry_error()) << err_description(__FILE__) << code_line(__LINE__);

            for (list<CGAL::Partition_traits_2<Kernel>::Polygon_2>::iterator it = partition_polys.begin();
                 it != partition_polys.end(); ++it) {
                obs.push_back(SimplePolygon(it->vertices_begin(), it->vertices_end()));
            }

            pts.clear();
        } else {
            Point p;
            std::stringstream strstr(line);
            strstr >> p;
            std::cout << p << std::endl;

            if (pts.empty() || p != pts.back())
                pts.push_back(p);
        }
    }

    f.close();

    Model model(0);
    model.set_obstacles(obs);
    model.write_obstacles_to_file("../data/utrecht_parsed.txt");
 */
}
