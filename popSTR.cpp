#include <iostream>
#include <string>

namespace computeReadAttributes    { int main(int argc, char const ** argv); }
namespace computePnSlippage        { int main(int argc, char const ** argv); }
namespace computePnSlippageDefault { int main(int argc, char const ** argv); }
namespace msGenotyper              { int main(int argc, char const ** argv); }
namespace msGenotyperDefault       { int main(int argc, char const ** argv); }


int main(int argc, char const ** argv)
{
  if (argc == 1)
  {
    std::cerr << "Usage: " << argv[0] << " <subcommand>\n\n"
              << "Where subcommand is one of: computeReadAttributes, computePnSlippage, computePnSlippageDefault, msGenotyper and msGenotyperDefault\n";
    return 0;
  }

  std::string subcommand = argv[1];

  if (subcommand == "computeReadAttributes")
    return computeReadAttributes::main(argc - 1, argv + 1);
  else if (subcommand == "computePnSlippage")
    return computePnSlippage::main(argc - 1, argv + 1);
  else if (subcommand == "computePnSlippageDefault")
    return computePnSlippageDefault::main(argc - 1, argv + 1);
  else if (subcommand == "msGenotyper")
    return msGenotyper::main(argc - 1, argv + 1);
  else if (subcommand == "msGenotyperDefault")
    return msGenotyperDefault::main(argc - 1, argv + 1);
  else
  {
    std::cerr << "Unrecognized subcommand, options are: "
              << "computeReadAttributes, computePnSlippage, computePnSlippageDefault, msGenotyper and "
              << "msGenotyperDefault.\n";
    return 1;
  }
}
