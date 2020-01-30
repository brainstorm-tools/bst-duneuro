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

#include <duneuro/common/matrix_utilities.hh>
//#include <duneuro/eeg/eeg_analytic_solution.hh>
#include <duneuro/io/dipole_reader.hh>
#include <duneuro/io/field_vector_reader.hh>
#include <duneuro/meeg/meeg_driver_factory.hh>

#include <iostream>

// compute the 2-norm of a vector
template <class T> T norm(const std::vector<T> &v) {
  return std::sqrt(std::inner_product(v.begin(), v.end(), v.begin(), T(0.0)));
}

template <class T>
T absolute_error(const std::vector<T> &num, const std::vector<T> &ana) {
  std::vector<T> diff;
  std::transform(num.begin(), num.end(), ana.begin(), std::back_inserter(diff),
                 [](const T &a, const T &b) { return a - b; });
  return norm(diff);
}

int run(const Dune::ParameterTree &config) {
  // set up driver
  auto driver = duneuro::MEEGDriverFactory<3>::make_meeg_driver(config);
  auto electrodes =
      duneuro::FieldVectorReader<double, 3>::read(config.sub("electrodes"));
  driver->setElectrodes(electrodes, config.sub("electrodes"));

  // read dipoles
  auto dipoles = duneuro::DipoleReader<double, 3>::read(config.sub("dipoles"));

  // store output in an output tree
  Dune::OutputTree output(config.get<std::string>("output.filename") + "." +
                          config.get<std::string>("output.extension"));

  // compute transfer matrix
  auto transfer = driver->computeEEGTransferMatrix(config.sub("solution"));

  auto solution = driver->makeDomainFunction();

  std::vector<double> absoluteErrors;

  // compute numerical solution transferred
  auto num_transfer =
      driver->applyEEGTransfer(*transfer, dipoles, config.sub("solution"));

  for (unsigned int i = 0; i < dipoles.size(); ++i) {

    // compute numerical solution direct
    driver->solveEEGForward(dipoles[i], *solution, config.sub("solution"));
    auto num_direct = driver->evaluateAtElectrodes(*solution);
    duneuro::subtract_mean(num_direct);

    // compute and store error measures
    auto prefix = std::string("dipole_") + std::to_string(i) + ".";
    absoluteErrors.push_back(absolute_error(num_direct, num_transfer[i]));
    output.set(prefix + "ae", absoluteErrors.back());
  }

  auto tolerance = config.get<double>("tolerance");
  for (unsigned int i = 0; i < absoluteErrors.size(); ++i) {
    if (absoluteErrors[i] > tolerance) {
      std::cerr << "error for source " << i
                << " is to big: " << absoluteErrors[i] << " > " << tolerance
                << std::endl;
      return -1;
    }
  }
  return 0;
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
    return run(config);
  } catch (Dune::Exception &e) {
    std::cerr << "Dune reported error: " << e << std::endl;
    return -1;
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
    return -1;
  }
}
