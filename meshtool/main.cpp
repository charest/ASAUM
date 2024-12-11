#include "driver.hpp"

#include "comm/mpi_comm.hpp"
#include "comm/comm_queue.hpp"
#include "mesh/mesh.hpp"
#include "utils/arguments.hpp"
#include "utils/mark.hpp"
#include "utils/thread.hpp"
#include "lua/lua_utils.hpp"

#include <cmath>
#include <string>

#ifdef ENABLE_FP_EXCEPTIONS
#include <fenv.h>
#endif

using namespace prl;

option_category_t OptionCat("GENERAL OPTIONS:");

static option_t<std::string> InputFile(
    "i",
    OptionCat,
    option_description_t("Input filename"),
    value_description_t("FILENAME"),
    required_t::yes);

// unused, just an example of a boolean input
static option_t<bool> Verbose(
    "v",
    OptionCat,
    option_description_t("Add verbosity"),
    required_t::no);

static option_t<std::string> LuaOption(
    "lua",
    OptionCat,
    option_description_t("Lua syntax to execute prior to loading the input file."));

#ifdef HAVE_THREADS
static option_t<bool> AsyncComm(
    "s",
    OptionCat,
    option_description_t("Use asyncronous communication."),
    required_t::no);
#endif

static option_t<bool> PrintTrace(
    "t",
    OptionCat,
    option_description_t("Print trace"),
    required_t::no);

static option_t<std::string> TraceFile(
    "trace-file",
    OptionCat,
    option_description_t("Dump trace to file"),
    value_description_t("FILENAME"));

////////////////////////////////////////////////////////////////////////////////
/// Main
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char ** argv) {
  
#ifdef ENABLE_FP_EXCEPTIONS
  feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
#endif


  // cant use parser yet, until mpi is initialized
  bool is_async = false;
  for (int i=0; i<argc; ++i)
    if (strcmp(argv[i], "-s") == 0)
      is_async = true;

  mpi_comm_t comm(argc, argv, is_async);
  
  MARK_FUNCTION();


  auto is_root = comm.is_root();
  
  // setup command line arguments
  auto ret = parse_command_line_arguments(argc, argv, is_root);

  if ( ret > 0 ) return consts::success;
  if ( ret < 0 ) return consts::failure;
  
  // create inputs type
  lua_t inputs(true, true);

  // set global variables
  if (LuaOption) {
    if (is_root) std::cout << "Executing command line lua..." << std::flush;
    inputs.run_string(LuaOption);
    if (is_root) std::cout << "done." << std::endl;
  }

  // get inputs
  std::string input_file_name = InputFile;
  
  // load input file
  inputs.loadfile( input_file_name );

  // set asyncronous comm
#ifdef HAVE_THREADS
  if (AsyncComm)
    comm_queue_t::set_async();
#endif
  
  // run the example
  driver_t drv(inputs, comm);
  drv.run();
  
  if (PrintTrace || TraceFile)
    DUMP_TRACE(comm, static_cast<bool>(PrintTrace), TraceFile);

  return consts::success;

}

