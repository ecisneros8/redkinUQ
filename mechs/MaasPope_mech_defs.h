#ifndef MECH_DEFS
#define MECH_DEFS

namespace mech
{

  //! Name
  const string mech_name = "MAASPOPE";

  //! Number of Reactions
  const int m_ii = 19;
  //! Number of Species
  const int m_kk = 9;
  //! Number of Elements
  const int m_mm = 3;

  //! Universal Gas Constant. [J/kmol/K]
  const double GasConstant = 8314.4621;

  //! One atmosphere [Pa]
  const double OneAtm = 1.01325e5;

  //! 1/3
  const double OneThird = 1.0/3.0;
  //! 5/16
  const double FiveSixteenths = 5.0/16.0;
  //! 1/6
  const double OneSixth = 1.0/6.0;
  //! sqrt(2)
  const double SqrtTwo = std::sqrt(2.0);

  

}

#endif
