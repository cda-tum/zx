#include "Simplify.hpp"

#include "Definitions.hpp"
#include "Rules.hpp"

#include <cstddef>
#include <iostream>

namespace zx {
    std::size_t simplifyVertices(ZXDiagram& diag, VertexCheckFun check,
                                 VertexRuleFun rule, const TerminationFun& isDone) {
        std::size_t n_simplifications = 0;
        bool        new_matches       = true;

        while (!isDone() && new_matches) {
            new_matches = false;
            for (auto [v, _]: diag.getVertices()) {
                if (check(diag, v)) {
                    rule(diag, v);
                    new_matches = true;
                    n_simplifications++;
                }
            }
        }

        return n_simplifications;
    }

    std::size_t simplifyEdges(ZXDiagram& diag, EdgeCheckFun check, EdgeRuleFun rule, const TerminationFun& isDone) {
        std::size_t n_simplifications = 0;
        bool        new_matches       = true;

        while (!isDone() && new_matches) {
            new_matches = false;
            for (auto [v0, v1]: diag.getEdges()) {
                if (diag.isDeleted(v0) || diag.isDeleted(v1) || !check(diag, v0, v1)) {
                    continue;
                }
                rule(diag, v0, v1);
                new_matches = true;
                n_simplifications++;
            }
        }

        return n_simplifications;
    }

    std::size_t gadgetSimp(ZXDiagram& diag, const TerminationFun& isDone) {
        std::size_t n_simplifications = 0;
        bool        new_matches       = true;

        while (!isDone() && new_matches) {
            new_matches = false;
            for (auto [v, _]: diag.getVertices()) {
                if (diag.isDeleted(v))
                    continue;

                if (checkAndFuseGadget(diag, v)) {
                    new_matches = true;
                    n_simplifications++;
                }
            }
        }
        return n_simplifications;
    }

    std::size_t idSimp(ZXDiagram& diag, const TerminationFun& isDone) {
        return simplifyVertices(diag, checkIdSimp, removeId, isDone);
    }

    std::size_t spiderSimp(ZXDiagram& diag, const TerminationFun& isDone) {
        return simplifyEdges(diag, checkSpiderFusion, fuseSpiders, isDone);
    }

    std::size_t localCompSimp(ZXDiagram& diag, const TerminationFun& isDone) {
        return simplifyVertices(diag, checkLocalComp, localComp, isDone);
    }

    std::size_t pivotPauliSimp(ZXDiagram& diag, const TerminationFun& isDone) {
        return simplifyEdges(diag, checkPivotPauli, pivotPauli, isDone);
    }

    std::size_t pivotSimp(ZXDiagram& diag, const TerminationFun& isDone) {
        return simplifyEdges(diag, checkPivot, pivot, isDone);
    }

    std::size_t interiorCliffordSimp(ZXDiagram& diag, const TerminationFun& isDone) {
        spiderSimp(diag, isDone);

        bool        new_matches       = true;
        std::size_t n_simplifications = 0;
        std::size_t n_id, n_spider, n_pivot, n_localComp;
        while (!isDone() && new_matches) {
            new_matches = false;
            n_id        = idSimp(diag, isDone);
            n_spider    = spiderSimp(diag, isDone);
            n_pivot     = pivotPauliSimp(diag, isDone);
            n_localComp = localCompSimp(diag, isDone);

            if (n_id + n_spider + n_pivot + n_localComp != 0) {
                new_matches = true;
                n_simplifications++;
            }
            // std::cout << "ID " << n_id << "\n";
            // std::cout << "SPIDER " << n_spider << "\n";
            // std::cout << "PIVOT PAULI" << n_pivot << "\n";
            // std::cout << "LOCALCOMP " << n_localComp << "\n";
        }
        return n_simplifications;
    }

    std::size_t cliffordSimp(ZXDiagram& diag, const TerminationFun& isDone) {
        bool        new_matches       = true;
        std::size_t n_simplifications = 0;
        std::size_t n_clifford, n_pivot;
        while (!isDone() && new_matches) {
            new_matches = false;
            n_clifford  = interiorCliffordSimp(diag, isDone);
            n_pivot     = pivotSimp(diag, isDone);
            if (n_clifford + n_pivot != 0) {
                new_matches = true;
                n_simplifications++;
            }
        }
        return n_simplifications;
    }

    std::size_t pivotgadgetSimp(ZXDiagram& diag, const TerminationFun& isDone) {
        return simplifyEdges(diag, checkPivotGadget, pivotGadget, isDone);
    }

    std::size_t fullReduce(ZXDiagram& diag, const TerminationFun& isDone) {
        diag.toGraphlike();
        interiorCliffordSimp(diag, isDone);

        // pivotgadgetSimp(diag);

        std::size_t n_gadget, n_pivot;
        std::size_t n_simplifications = 0;
        while (!isDone()) {
            cliffordSimp(diag, isDone);
            n_gadget = gadgetSimp(diag, isDone);
            interiorCliffordSimp(diag, isDone);
            n_pivot = pivotgadgetSimp(diag, isDone);
            if (n_gadget + n_pivot == 0)
                break;
            n_simplifications += n_gadget + n_pivot;
        }
        diag.removeDisconnectedSpiders();

        return n_simplifications;
    }

    std::size_t fullReduceApproximate(ZXDiagram& diag, fp tolerance, const TerminationFun& isDone) {
        auto        nSimplifications = fullReduce(diag, isDone);
        std::size_t newSimps         = 0;
        do {
            diag.approximateCliffords(tolerance);
            newSimps = fullReduce(diag, isDone);
            nSimplifications += newSimps;
        } while (!isDone() && newSimps);
        return nSimplifications;
    }
} // namespace zx
