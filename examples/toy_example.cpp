#include <sbg/interval.hpp>
#include <sbg/pw_map.hpp>
#include <sbg/sbg.hpp>


void lets_define_a_sbg_graph()
{
  SBG::LIB::CanonSBG graph;
  SBG::LIB::OrdSet s1;
  // for (int j = 1; j <= 20; j++) {
  //   SBG::LIB::Interval i(j*100+1, 1, (j+1)*100);
  //   s1.emplace_hint(s1.end(), i);
  // }
  SBG::LIB::Interval i_1(1, 1, 5);
  SBG::LIB::Interval i_2(6, 1, 10);
  SBG::LIB::Interval i_3(11, 1, 15);
  SBG::LIB::Interval i_4(16, 1, 20);
  s1.emplace_hint(s1.end(), i_1);
  s1.emplace_hint(s1.end(), i_2);
  s1.emplace_hint(s1.end(), i_3);
  s1.emplace_hint(s1.end(), i_4);
  graph.set_V(s1);

  SBG::LIB::CanonPWMap ms1;

  SBG::LIB::LExp le0(0, 20);
  ms1.emplace(SBG::LIB::CanonMap(SBG::LIB::Interval(1, 1, 1), le0));

  SBG::LIB::LExp le1(1, 10);
  ms1.emplace(SBG::LIB::CanonMap(i_1, le1));

  SBG::LIB::LExp le2(1, 9);
  ms1.emplace(SBG::LIB::CanonMap(SBG::LIB::Interval(2, 1, 5), le2));

  SBG::LIB::LExp le3(0, 15);
  ms1.emplace(SBG::LIB::CanonMap(SBG::LIB::Interval(6, 1, 10), le3));

  ms1.emplace(SBG::LIB::CanonMap(SBG::LIB::Interval(7, 1, 10), le2));

  graph.set_map1(ms1);

  SBG::LIB::CanonPWMap ms2;

  SBG::LIB::Interval int_1(11, 1, 15);
  SBG::LIB::LExp le4 = SBG::LIB::LExp(1, 0) - SBG::LIB::LExp(0, 10);
  ms2.emplace(SBG::LIB::CanonMap(int_1, le4));
  SBG::LIB::LExp le5 = SBG::LIB::LExp(1, 0) - SBG::LIB::LExp(0, 9);
  ms2.emplace(SBG::LIB::CanonMap(int_1, le5));

  SBG::LIB::Interval int_2(15, 1, 15);
  SBG::LIB::LExp le5_6_1(0, 7);
  SBG::LIB::LExp le5_6_2(0, 8);
  SBG::LIB::LExp le5_6_3(0, 9);
  SBG::LIB::LExp le5_6_4(0, 10);
  ms2.emplace(SBG::LIB::CanonMap(int_2,  le5_6_1));
  ms2.emplace(SBG::LIB::CanonMap(int_2,  le5_6_2));
  ms2.emplace(SBG::LIB::CanonMap(int_2,  le5_6_3));
  ms2.emplace(SBG::LIB::CanonMap(int_2,  le5_6_4));

  SBG::LIB::Interval int_3(16, 1, 19);
  SBG::LIB::LExp le6 = SBG::LIB::LExp(1, 0) - SBG::LIB::LExp(0, 9);
  ms2.emplace(SBG::LIB::CanonMap(int_2, le6));

  SBG::LIB::Interval int_4(20, 1, 20);
  SBG::LIB::LExp le7 = SBG::LIB::LExp(0, 1);
  ms2.emplace(SBG::LIB::CanonMap(int_4, le7));

  graph.set_map2(ms2);

  std::cout << graph << std::endl;
}
