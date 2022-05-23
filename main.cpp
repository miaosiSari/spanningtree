#include "vertex.h"

int main(void){
    helper h;
    std::unordered_map<size_t, std::shared_ptr<vertex>> adj;
    Eigen::MatrixXf m(3, 3);
    m(0, 0) = m(1, 1) = m(2, 2) = 2;
    m(0, 1) = m(1, 0) = m(0, 2) = m(2, 0) = m(1, 2) = m(2, 1) = -1;
    auto r = h.construct(m, 3, adj);
    h.print(adj);
    std::cout << h.Kirchhoff(r) << " " << h.Kirchhoff(m) << "\n";
    std::cout << *h.genKn(6) << " " << h.Kirchhoff(h.genKn(6)) << "\n";
    std::cout << *h.construct(h.genKn(6), 3, adj) << " " <<  (int)h.Kirchhoff(h.construct(h.genKn(6), 3, adj)) << "\n";
    std::cout << *h.genKmn(3, 3) << " " << h.Kirchhoff(h.genKmn(3, 3)) << "\n";
    std::cout << *h.construct(h.genKmn(3, 3), 3, adj) << " " <<  (int)h.Kirchhoff(h.construct(h.genKmn(3, 3), 3, adj)) << "\n";
    auto rand = *h.genrandom(5);
    std::cout << rand << " " << h.Kirchhoff(rand) << " " << h.Kirchhoff(h.construct(rand, 3, adj)) << "\n";
}