/*
 * This file is part of JKQ QCEC library which is released under the MIT license.
 * See file README.md or go to http://iic.jku.at/eda/research/quantum_verification/ for more information.
 */

#include <ImprovedDDEquivalenceChecker.hpp>

namespace ec {
    qc::MatrixDD ImprovedDDEquivalenceChecker::createInitialMatrix() {
        auto e = dd->makeIdent(nqubits);
        dd->incRef(e);

        std::vector<bool> ancillary(nqubits);
        for (auto q = static_cast<dd::Qubit>(nqubits - 1); q >= 0; --q) {
            if (qc1.logicalQubitIsAncillary(q) && qc2.logicalQubitIsAncillary(q)) {
                bool found1  = false;
                bool isidle1 = false;
                for (const auto& in1: initial1) {
                    if (in1.second == q) {
                        found1  = true;
                        isidle1 = qc1.isIdleQubit(in1.first);
                        break;
                    }
                }
                bool found2  = false;
                bool isidle2 = false;
                for (const auto& in2: initial2) {
                    if (in2.second == q) {
                        found2  = true;
                        isidle2 = qc2.isIdleQubit(in2.first);
                        break;
                    }
                }

                // qubit only really exists or is acted on in one of the circuits
                if ((found1 ^ found2) || (isidle1 ^ isidle2)) {
                    ancillary[q] = true;
                }
            }
        }
        e = dd->reduceAncillae(e, ancillary);
        return e;
    }

    qc::MatrixDD ImprovedDDEquivalenceChecker::createGoalMatrix() {
        auto goalMatrix = dd->makeIdent(nqubits);
        dd->incRef(goalMatrix);
        goalMatrix = dd->reduceAncillae(goalMatrix, ancillary2, RIGHT);
        goalMatrix = dd->reduceGarbage(goalMatrix, garbage2, RIGHT);
        goalMatrix = dd->reduceAncillae(goalMatrix, ancillary1, LEFT);
        goalMatrix = dd->reduceGarbage(goalMatrix, garbage1, LEFT);
        return goalMatrix;
    }



    dd::fp* ImprovedDDEquivalenceChecker::traceRecur(dd::Package::mNode * cur, std::map<dd::Package::mNode *, dd::fp *> *Node_Table) {
        // std::cout<<0<<std::endl;
        std::map<dd::Package::mNode *, dd::fp *>::iterator it;
        if (!(dd::Package::mNode::isTerminal(cur->e[0].p)) && (Node_Table->find(cur->e[0].p) == Node_Table->end()))
        {
            // std::cout<<1<<std::endl;
            dd::fp * dum = traceRecur(cur->e[0].p, Node_Table);
            // std::cout<<2<<std::endl;
        }
        if (!(dd::Package::mNode::isTerminal(cur->e[3].p)) && (Node_Table->find(cur->e[3].p) == Node_Table->end()))
            dd::fp * dum = traceRecur(cur->e[3].p, Node_Table);

        dd::fp * left;
        dd::fp * right;
        if (dd::Package::mNode::isTerminal(cur->e[0].p)) {
            left = new dd::fp[2]; 
            left[0] = 1.0;
            left[1] = 0.0;
        }
        else {
            it = Node_Table->find(cur->e[0].p);
            assert(it != Node_Table->end());
            left = it->second;
            // std::cout << "  left: " << left << std::endl;
        }
        if (dd::Package::mNode::isTerminal(cur->e[3].p)) {
            right = new dd::fp[2]; 
            right[0] = 1.0;
            right[1] = 0.0;
        }
        else {
            it = Node_Table->find(cur->e[3].p);
            assert(it != Node_Table->end());
            right = it->second;
            // std::cout << "  right: " << right << std::endl;
        }
        // std::cout << "  wl: " << cur->e[0].w << std::endl;
        // std::cout << "  wr: " << cur->e[3].w << std::endl;
        dd::fp * w0 = complex2fp(cur->e[0].w);
        dd::fp * w1 = complex2fp(cur->e[3].w);
        dd::fp * le = mul(w0, left);
        dd::fp * re = mul(w1, right);
        // std::cout << "  le: " << le << std::endl;
        // std::cout << "  re: " << re << std::endl;
        dd::fp * sum = add(le, re);
        (*Node_Table)[cur] = sum;
        // std::cout << "  sum: " << sum << std::endl << std::endl;
        // clear 
        if (dd::Package::mNode::isTerminal(cur->e[0].p)) delete [] left;
        if (dd::Package::mNode::isTerminal(cur->e[3].p)) delete [] right;
        delete [] w0;
        delete [] w1;
        delete [] le;
        delete [] re;
        
        return sum;
    }

    dd::fp * ImprovedDDEquivalenceChecker::complex2fp(const dd::Complex &c)
    {
        dd::fp * a = new dd::fp[2];
        a[0] = dd::CTEntry::val(c.r);
        a[1] = dd::CTEntry::val(c.i);
        return a;
    }

    dd::fp * ImprovedDDEquivalenceChecker::add(dd::fp * a, dd::fp * b)
    {
        dd::fp * r = new dd::fp[2];
        r[0] = a[0] + b[0];
        r[1] = a[1] + b[1];
        return r;
    }

    dd::fp * ImprovedDDEquivalenceChecker::mul(dd::fp * a, dd::fp * b)
    {
        dd::fp * r = new dd::fp[2];
        r[0] = a[0] * b[0] - a[1] * b[1];
        r[1] = a[0] * b[1] + a[1] * b[0];
        return r;
    }

    /// Use dedicated method to check the equivalence of both provided circuits
    EquivalenceCheckingResults ImprovedDDEquivalenceChecker::check(const Configuration& config) {
        EquivalenceCheckingResults results{};
        setupResults(results);
        results.strategy = config.strategy;

        auto start = std::chrono::steady_clock::now();
        runPreCheckPasses(config);
        auto endPreprocessing = std::chrono::steady_clock::now();

        auto perm1     = initial1;
        auto perm2     = initial2;
        results.result = createInitialMatrix();

        switch (config.strategy) {
            case ec::Strategy::Naive:
                checkNaive(results.result, perm1, perm2);
                break;
            case ec::Strategy::Proportional:
                checkProportional(results.result, perm1, perm2);
                break;
            case ec::Strategy::Lookahead:
                checkLookahead(results.result, perm1, perm2);
                break;
            default:
                throw std::invalid_argument("Strategy " + toString(config.strategy) + " not supported by ImprovedDDEquivalenceChecker");
        }

        // finish first circuit
        while (it1 != end1) {
            applyGate(qc1, it1, results.result, perm1, LEFT);
            ++it1;
        }

        //finish second circuit
        while (it2 != end2) {
            applyGate(qc2, it2, results.result, perm2, RIGHT);
            ++it2;
        }

        qc::QuantumComputation::changePermutation(results.result, perm1, output1, dd, LEFT);
        qc::QuantumComputation::changePermutation(results.result, perm2, output2, dd, RIGHT);
        results.result = dd->reduceGarbage(results.result, garbage1, LEFT);
        results.result = dd->reduceGarbage(results.result, garbage2, RIGHT);
        results.result = dd->reduceAncillae(results.result, ancillary1, LEFT);
        results.result = dd->reduceAncillae(results.result, ancillary2, RIGHT);

        // fidelity
        if (isFid)
        {
            std::map<dd::Package::mNode *, dd::fp *> Node_Table;
            std::map<dd::Package::mNode *, dd::fp *>::iterator it;
            dd::fp * sum = traceRecur(results.result.p, &Node_Table);
            dd::fp * w = complex2fp(results.result.w);
            dd::fp * sum_w = mul(w, sum);
            // std::cout << "  sum_w: " << sum_w << std::endl;
            fid = (sum_w[0]*sum_w[0] + sum_w[1]*sum_w[1])/(pow(2, 2*nqubits));
            std::cout << "  Fidelity: " << fid << std::endl;
            // clear
            delete [] w;
            delete [] sum_w;
            for (it = Node_Table.begin(); it != Node_Table.end(); it++) delete [] it->second;
            Node_Table.clear();
        }

        results.equivalence = equals(results.result, createGoalMatrix());
        results.maxActive   = std::max(results.maxActive, dd->mUniqueTable.getMaxActiveNodes());

        auto                          endVerification   = std::chrono::steady_clock::now();
        std::chrono::duration<double> preprocessingTime = endPreprocessing - start;
        std::chrono::duration<double> verificationTime  = endVerification - endPreprocessing;
        results.preprocessingTime                       = preprocessingTime.count();
        results.verificationTime                        = verificationTime.count();

        return results;
    }

    /// Alternate between LEFT and RIGHT applications
    void ImprovedDDEquivalenceChecker::checkNaive(qc::MatrixDD& result, qc::Permutation& perm1, qc::Permutation& perm2) {
        while (it1 != end1 && it2 != end2) {
            applyGate(qc1, it1, result, perm1, LEFT);
            ++it1;
            applyGate(qc2, it2, result, perm2, RIGHT);
            ++it2;
        }
    }

    /// Alternate according to the gate count ratio between LEFT and RIGHT applications
    void ImprovedDDEquivalenceChecker::checkProportional(qc::MatrixDD& result, qc::Permutation& perm1, qc::Permutation& perm2) {
        auto ratio  = static_cast<unsigned int>(std::round(
                static_cast<double>(std::max(qc1.getNops(), qc2.getNops())) /
                static_cast<double>(std::min(qc1.getNops(), qc2.getNops()))));
        auto ratio1 = (qc1.getNops() > qc2.getNops()) ? ratio : 1;
        auto ratio2 = (qc1.getNops() > qc2.getNops()) ? 1 : ratio;

        while (it1 != end1 && it2 != end2) {
            for (unsigned int i = 0; i < ratio1 && it1 != end1; ++i) {
                applyGate(qc1, it1, result, perm1, LEFT);
                ++it1;
            }
            for (unsigned int i = 0; i < ratio2 && it2 != end2; ++i) {
                applyGate(qc2, it2, result, perm2, RIGHT);
                ++it2;
            }
        }
    }

    /// Look-ahead LEFT and RIGHT and choose the more promising option
    void ImprovedDDEquivalenceChecker::checkLookahead(qc::MatrixDD& result, qc::Permutation& perm1, qc::Permutation& perm2) {
        qc::MatrixDD left{}, right{}, saved{};
        bool         cachedLeft = false, cachedRight = false;

        while (it1 != end1 && it2 != end2) {
            if (!cachedLeft) {
                // stop if measurement is encountered
                if ((*it1)->getType() == qc::Measure)
                    break;

                auto nq = (*it1)->getNqubits();
                (*it1)->setNqubits(nqubits);
                left = (*it1)->getDD(dd, perm1);
                dd->incRef(left);
                (*it1)->setNqubits(nq);
                ++it1;
                cachedLeft = true;
            }

            if (!cachedRight) {
                // stop if measurement is encountered
                if ((*it2)->getType() == qc::Measure)
                    break;

                auto nq = (*it2)->getNqubits();
                (*it2)->setNqubits(nqubits);
                right = (*it2)->getInverseDD(dd, perm2);
                dd->incRef(right);
                (*it2)->setNqubits(nq);
                ++it2;
                cachedRight = true;
            }

            saved          = result;
            auto lookLeft  = dd->multiply(left, saved);
            auto lookRight = dd->multiply(saved, right);

            auto nc1 = dd->size(lookLeft);
            auto nc2 = dd->size(lookRight);

            if (nc1 <= nc2) {
                result = lookLeft;
                dd->decRef(left);
                cachedLeft = false;
            } else {
                result = lookRight;
                dd->decRef(right);
                cachedRight = false;
            }
            dd->incRef(result);
            dd->decRef(saved);
            dd->garbageCollect();
        }

        if (cachedLeft) {
            saved  = result;
            result = dd->multiply(left, saved);
            dd->incRef(result);
            dd->decRef(saved);
            dd->decRef(left);
            dd->garbageCollect();
        }

        if (cachedRight) {
            saved  = result;
            result = dd->multiply(saved, right);
            dd->incRef(result);
            dd->decRef(saved);
            dd->decRef(right);
            dd->garbageCollect();
        }
    }
} // namespace ec
