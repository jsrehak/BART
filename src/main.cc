#include <memory>
#include <iostream>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/mpi.h>

#include "problem/parameters_dealii_handler.h"
#include "problem/locator.h"

int main(int argc, char* argv[]) {
  try {
    if (argc != 2) {
      std::cerr
          << "Call the program as mpirun -np num_proc xtrans input_file_name"
          << std::endl;
      return 1;
    }
    dealii::ParameterHandler prm;

    // New parameters handler
    auto parameter_handler_ptr =
        std::make_shared<bart::problem::ParametersDealiiHandler>();

    // Declare input strings, declare both using the new handler, and the old
    // method, as not all have been moved (and may not be moved)
    parameter_handler_ptr->SetUp(prm);

    const std::string filename{argv[1]};
    prm.parse_input(filename, "");

    parameter_handler_ptr->Parse(prm);
    bart::problem::Locator::Provide(parameter_handler_ptr);

    dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  } catch (std::exception &exc) {
    std::cerr << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  } catch (...) {
    std::cerr << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }
  return 0;
}
