#pragma once

#include "Rules.hpp"
#include "ZXDiagram.hpp"

namespace zx {
    using VertexCheckFun                          = decltype(checkIdSimp);
    using VertexRuleFun                           = decltype(removeId);
    using EdgeCheckFun                            = decltype(checkSpiderFusion);
    using EdgeRuleFun                             = decltype(fuseSpiders);
    const std::function<bool(void)> defaultIsDone = []() { return false; };
    using TerminationFun                          = decltype(defaultIsDone);

    std::size_t simplifyVertices(ZXDiagram& diag, VertexCheckFun check,
                                 VertexRuleFun rule, const TerminationFun& isDone = defaultIsDone);

    std::size_t simplifyEdges(ZXDiagram& diag, EdgeCheckFun check,
                              EdgeRuleFun rule, const TerminationFun& isDone = defaultIsDone);

    std::size_t gadgetSimp(ZXDiagram& diag, const TerminationFun& isDone = defaultIsDone);

    std::size_t idSimp(ZXDiagram& diag, const TerminationFun& isDone = defaultIsDone);

    std::size_t spiderSimp(ZXDiagram& diag, const TerminationFun& isDone = defaultIsDone);

    std::size_t localCompSimp(ZXDiagram& diag, const TerminationFun& isDone = defaultIsDone);

    std::size_t pivotPauliSimp(ZXDiagram& diag, const TerminationFun& isDone = defaultIsDone);

    std::size_t pivotSimp(ZXDiagram& diag, const TerminationFun& isDone = defaultIsDone);

    std::size_t interiorCliffordSimp(ZXDiagram& diag, const TerminationFun& isDone = defaultIsDone);

    std::size_t cliffordSimp(ZXDiagram& diag, const TerminationFun& isDone = defaultIsDone);

    std::size_t pivotgadgetSimp(ZXDiagram& diag, const TerminationFun& isDone = defaultIsDone);

    std::size_t fullReduce(ZXDiagram& diag, const TerminationFun& isDone = defaultIsDone);
    std::size_t fullReduceApproximate(ZXDiagram& diag, fp tolerance, const TerminationFun& isDone = defaultIsDone);

} // namespace zx
