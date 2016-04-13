#include <libplugin/plugin.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include "rhf/rhf.hpp"

INIT_PLUGIN

extern "C" 
int read_options(std::string name, psi::Options& options) { return false; }


/* PLUGIN MAIN() */
extern "C" 
psi::PsiReturnType plugin_main(psi::Options& options)
{

  /* Your code goes here */
  /* BEGIN SAMPLE */
  boost::shared_ptr<psi::Molecule> mol = psi::Process::environment.molecule();
  rhf::RHF rhf(mol, options);
  double E = rhf.compute_energy();

  /* END SAMPLE */

  return psi::Success;
}
