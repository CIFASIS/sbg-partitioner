#include <fstream>
#include <iostream>

namespace sbg_partitioner {

/**
 * @brief Write the input in a given file, useful for regression tests.
 * 
 * @tparam T Input type
 * @param filename output file
 * @param input object to be logged
 */
template<typename T>
void logger(const std::string& filename, const T& input)
{
  // For testing purposes
  std::ofstream out(filename);
  std::streambuf *coutbuf = std::cout.rdbuf();
  std::cout.rdbuf(out.rdbuf());
  std::cout << input;
  std::cout.rdbuf(coutbuf);
}
}