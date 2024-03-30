#include "generators/boost_random_number_generator.hpp"
#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"
#include "generators/known_polytope_generators.h"
#include "generators/h_polytopes_generator.h"
#include "random_walks/random_walks.hpp"
#include "random_walks/gaussian_hamiltonian_monte_carlo_exact_walk.hpp"
#include "random_walks/gaussian_rdhr_walk.hpp"
#include "random_walks/hamiltonian_monte_carlo_walk.hpp"
#include "volume/volume_cooling_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"

typedef Cartesian<double> Kernel;
typedef typename Kernel::Point Point;
typedef BoostRandomNumberGenerator<boost::random::mt19937, double> RNGType;
typedef HPolytope<Point> HPolytopeType;

int main() {
	// Generating a random h polytope of dim =3 with 6 segments.
	HPolytopeType HP = random_hpoly<HPolytopeType, boost::random::mt19937>(3, 6);
	std::cout<<"Polytope: \n";
	HP.print();
	std::cout<<"\n";

	// Setup parameters for calculating volume
	int walk_len = 10 + HP.dimension()/10;
	double e = 0.05;

	// Calculating volume of the passed polytope
	double volume = volume_cooling_balls
        <BallWalk, RNGType, HPolytopeType>(HP, e, walk_len).second;

	// Since CG works well on Hpoly with dim > 200
	HPolytopeType cube = generate_cube<HPolytopeType>(210,false);
	walk_len = 10 + cube.dimension()/10;
	// Cooling gaussians with rdhr
	double volume_cooling_gaussians_rdhr = volume_cooling_gaussians<GaussianHamiltonianMonteCarloExactWalk, RNGType, HPolytopeType>(cube,e, walk_len);
	// cooling gaussians with hmc
	// double volume_cooling_gaussians_hmc = volume_cooling_gaussians<>();
	// // cooling gaussians with exact hmc
	// double volume_cooling_gaussians_exact_hmc = volume_cooling_gaussians<>();
    std::cout << "Volume of the h-polytope: " << volume << std::endl;
	std::cout << "Volume of the random h-polytope using cg with HMC: " << volume_cooling_gaussians_rdhr << std::endl;

	return 0;
}
