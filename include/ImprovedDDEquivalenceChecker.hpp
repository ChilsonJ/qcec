/*
 * This file is part of JKQ QCEC library which is released under the MIT license.
 * See file README.md or go to http://iic.jku.at/eda/research/quantum_verification/ for more information.
 */

#ifndef QUANTUMCIRCUITEQUIVALENCECHECKING_IMPROVEDDDEQUIVALENCECHECKER_HPP
#define QUANTUMCIRCUITEQUIVALENCECHECKING_IMPROVEDDDEQUIVALENCECHECKER_HPP

#include <memory>
#include <array>
#include <unordered_set>
#include <chrono>

#include "EquivalenceChecker.hpp"

namespace ec {

	class ImprovedDDEquivalenceChecker: public EquivalenceChecker {

		/// Alternate between LEFT and RIGHT applications
		void checkNaive(dd::Edge& result, qc::permutationMap& perm1, qc::permutationMap& perm2);
		/// Alternate according to the gate count ratio between LEFT and RIGHT applications
		void checkProportional(dd::Edge& result, qc::permutationMap& perm1, qc::permutationMap& perm2);
		/// Look-ahead LEFT and RIGHT and choose the more promising option
		void checkLookahead(dd::Edge& result, qc::permutationMap& perm1, qc::permutationMap& perm2);

	protected:
		/// Create the initial matrix used for the G->I<-G' scheme.
		/// [1 0] if the qubit is no ancillary or it is acted upon by both circuits
		/// [0 1]
		///
		/// [1 0] for an ancillary that is present in one circuit and not acted upon in the other
		/// [0 0]
		/// \return initial matrix
		dd::Edge createInitialMatrix();

		/// Create the goal matrix used for the G->I<-G' scheme.
		/// [1 0] if the qubit is no ancillary
		/// [0 1]
		///
		/// [1 0] for an ancillary that is present in either circuit
		/// [0 0]
		/// \return goal matrix
		dd::Edge createGoalMatrix();

	public:
		ImprovedDDEquivalenceChecker(qc::QuantumComputation& qc1, qc::QuantumComputation& qc2):
				EquivalenceChecker(qc1, qc2) {
			method = Method::G_I_Gp;
		}

		/// Use dedicated method to check the equivalence of both provided circuits
		EquivalenceCheckingResults check(const Configuration& config) override;
		EquivalenceCheckingResults check() override { return check(Configuration{}); }
	};

}

#endif //QUANTUMCIRCUITEQUIVALENCECHECKING_IMPROVEDDDEQUIVALENCECHECKER_HPP
