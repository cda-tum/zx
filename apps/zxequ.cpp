#include "Definitions.hpp"
#include "QuantumComputation.hpp"
#include "Rules.hpp"
#include "Simplify.hpp"
#include "ZXDiagram.hpp"

#include <algorithm>
#include <ctime>
#include <iostream>

void print_diag(zx::ZXDiagram d0) {
      std::cout
      << d0.getNEdges() << "\n";
  std::cout << d0.getNVertices() << "\n";

  for (auto [from, to] : d0.getEdges()) {
    std::cout << from
              << (d0.getEdge(from, to).value().type ==
              zx::EdgeType::Hadamard
                      ? "- -"
                      : "---")
              << to << "\n";
  }
  std::cout << ""
            << "\n";

  for (int i = 0; i < d0.getInputs().size(); i++) {
    std::cout << d0.getInputs()[i] << "--" << d0.getOutputs()[i] << "\n";
  }
  std::cout << ""
            << "\n";

  for (auto [v, data] : d0.getVertices())
    std::cout << v << " p: " << data.phase <<", q:" << ((int)data.qubit) <<
    ", r:" << (data.col)<<"\n";
  std::cout << ""
            << "\n";
  for (auto [v, data] : d0.getVertices()) {
    std::cout << v << " p:" << data.phase << " boundary "
              << (data.type == zx::VertexType::Boundary ? "True" : "False")
              << " type " << (d0.type(v) == zx::VertexType::Z ? "Z" : "X")
              << "\n";
  }
  
}
int main(int argc, char **argv) {
  qc::QuantumComputation c0(argv[1]);
  qc::QuantumComputation c1(argv[2]);
  // std::cout << "ANCILLAE: " << ((int)c1.getNancillae())<<"\n";
  zx::ZXDiagram d0(c0);
  zx::ZXDiagram d1(c1);
  d0.invert();
  d0.concat(d1);
  // zx::cliffordSimp(d0);
  // zx::spiderSimp(d0);

  std::clock_t c_start = std::clock();
  zx::fullReduce(d0);
  qc::Permutation& p0 = c0.outputPermutation;
  qc::Permutation& p1 = c1.outputPermutation;
  qc::Permutation p;
  for(auto& [from_0, to_0]: p0){
    auto k_v = std::find_if(p1.begin(), p1.end(), [&](auto& k_v) {return k_v.second==to_0;});
    p[from_0] = k_v->first;
  }
  std::clock_t c_end = std::clock();

  double time_elapsed_s = (c_end-c_start) / (double CLOCKS_PER_SEC);
  // std::cout << "CPU time used: " << time_elapsed_ms << " ms\n";


  std::cout << static_cast<int>(c0.getNqubits()) << ", " << c0.getNops() <<
    ", " << c1.getNops() << ", " << time_elapsed_s << ", ";

  print_diag(d0);
  // std::cout << "" << "\n";

  // if(d0.isIdentity(p))
  //   std::cout << "IDENTITY " << "\n";

  if (d0.isIdentity(p)) {
    std::cout << "TRUE";
  } else {
    std::cout << "FALSE";
  }
  // if (d0.getInputs().size() == d0.getNEdges()) {
  //   std::cout << "TRUE";
  // } else {
  //   std::cout << "FALSE";
  // }
  return 0;
}