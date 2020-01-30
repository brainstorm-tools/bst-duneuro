// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <algorithm>

#include <dune/common/exceptions.hh>         // We use exceptions
#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>

#include <dune/testtools/outputtree.hh>
#include <src/include/simbiosphere/analytic_solution.hh>

#include <duneuro/io/dipole_reader.hh>
#include <duneuro/io/field_vector_reader.hh>
#include <duneuro/meeg/meeg_driver_factory.hh>

#include <iostream>

// compute the 2-norm of a vector
template <class T> T norm(const std::vector<T> &v) {
  return std::sqrt(std::inner_product(v.begin(), v.end(), v.begin(), T(0.0)));
}

// compute \|num-ana\|/\|ana\|
template <class T>
T relative_error(const std::vector<T> &num, const std::vector<T> &ana) {
  std::vector<T> diff;
  std::transform(num.begin(), num.end(), ana.begin(), std::back_inserter(diff),
                 [](const T &a, const T &b) { return a - b; });
  return norm(diff) / norm(ana);
}

// compute \|num\|/\|ana\|
template <class T>
T magnitude_error(const std::vector<T> &num, const std::vector<T> &ana) {
  return norm(num) / norm(ana);
}

// compute \| num/\|num\| - ana/\|ana\|
template <class T>
T rdm_error(const std::vector<T> &num, const std::vector<T> &ana) {
  auto nn = norm(num);
  auto na = norm(ana);
  std::vector<T> diff;
  std::transform(num.begin(), num.end(), ana.begin(), std::back_inserter(diff),
                 [nn, na](const T &a, const T &b) { return a / nn - b / na; });
  return norm(diff);
}

// subtract the mean of each entry
template <class T> void subtract_mean(std::vector<T> &sol) {
  T mean = std::accumulate(sol.begin(), sol.end(), T(0.0)) / sol.size();
  for (auto &s : sol)
    s -= mean;
}

void run(const Dune::ParameterTree &config) {
  // set up driver
  auto driver = duneuro::MEEGDriverFactory<3>::make_meeg_driver(config);
  auto electrodes =
      duneuro::FieldVectorReader<double, 3>::read(config.sub("electrodes"));
  driver->setElectrodes(electrodes, config.sub("electrodes"));

  // read dipoles
  auto dipoles = duneuro::DipoleReader<double, 3>::read(config.sub("dipoles"));

  // read gemetric information

  auto radii = config.get<std::array<double, 4>>("analytic_solution.radii");
  auto center = config.get<std::array<double, 3>>("analytic_solution.center");
  auto conductivities =
      config.get<std::array<double, 4>>("analytic_solution.conductivities");

  // create storage for solution
  auto solution = driver->makeDomainFunction();

  // store output in an output tree
  Dune::OutputTree output(config.get<std::string>("output.filename") + "." +
                          config.get<std::string>("output.extension"));

  for (unsigned int i = 0; i < dipoles.size(); ++i) {
    // compute numerical solution
    driver->solveEEGForward(dipoles[i], *solution, config.sub("solution"));
    auto num = driver->evaluateAtElectrodes(*solution);
    subtract_mean(num);

    auto dipole = dipoles[i];

    std::vector<std::array<double, 3>> elec;
    std::array<double, 3> tmp;
    for (auto i : electrodes) {
      std::copy_n(i.begin(), 3, tmp.begin());
      elec.push_back(tmp);
    }

    std::array<double, 3> dipposition;
    std::copy_n(dipole.position().begin(), 3, dipposition.begin());

    std::array<double, 3> dipmoment;
    std::copy_n(dipole.moment().begin(), 3, dipmoment.begin());

    // compute analytic solution
    auto ana = simbiosphere::analytic_solution(radii, center, conductivities,
                                               elec, dipposition, dipmoment);
    subtract_mean(ana);

    // compute and store error measures
    auto prefix = std::string("dipole_") + std::to_string(i) + ".";
    output.set(prefix + "re", relative_error(num, ana));
    output.set(prefix + "mag", magnitude_error(num, ana));
    output.set(prefix + "rdm", rdm_error(num, ana));
  }
}

int main(int argc, char **argv) {
  try {
    // Maybe initialize MPI
    Dune::MPIHelper::instance(argc, argv);
    if (argc != 2) {
      std::cerr << "please provide a config file";
      return -1;
    }
    Dune::ParameterTree config;
    Dune::ParameterTreeParser::readINITree(argv[1], config);
    run(config);
  } catch (Dune::Exception &e) {
    std::cerr << "Dune reported error: " << e << std::endl;
    return -1;
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
    return -1;
  }
}
