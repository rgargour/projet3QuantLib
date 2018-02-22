#pragma once

/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2002, 2003, 2004 Ferdinando Ametrano
Copyright (C) 2002, 2003 RiskMap srl
Copyright (C) 2003, 2004, 2005, 2007 StatPro Italia srl
Copyright (C) 2007 Affine Group Limited
This file is part of QuantLib, a free-software/open-source library
for financial quantitative analysts and developers - http://quantlib.org/
QuantLib is free software: you can redistribute it and/or modify it
under the terms of the QuantLib license.  You should have received a
copy of the license along with this program; if not, please email
<quantlib-dev@lists.sf.net>. The license is also available online at
<http://quantlib.org/license.shtml>.
This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

/*! \file binomialengine.hpp
\brief Binomial option engine
*/

#ifndef binomial_engine_hpp
#define binomial_engine_hpp

#include <ql/methods/lattices/binomialtree.hpp>
#include <ql/methods/lattices/bsmlattice.hpp>
#include <ql/math/distributions/normaldistribution.hpp>
#include <ql/pricingengines/vanilla/discretizedvanillaoption.hpp>
#include <ql/pricingengines/greeks.hpp>
#include <ql/processes/blackscholesprocess.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/termstructures/volatility/equityfx/blackconstantvol.hpp>

namespace QuantLib {

	//! Pricing engine for vanilla options using binomial trees
	/*! \ingroup vanillaengines
	\test the correctness of the returned values is tested by
	checking it against analytic results.
	\todo Greeks are not overly accurate. They could be improved
	by building a tree so that it has three points at the
	current time. The value would be fetched from the middle
	one, while the two side points would be used for
	estimating partial derivatives.
	*/
	template <class T>
	class BinomialVanillaEngine_2 : public VanillaOption::engine {
	public:
		BinomialVanillaEngine_2(
			const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
			Size timeSteps)
			: process_(process), timeSteps_(timeSteps) {
			QL_REQUIRE(timeSteps >= 2,
				"at least 2 time steps required, "
				<< timeSteps << " provided");
			registerWith(process_);
		}
		void calculate() const;
	private:
		boost::shared_ptr<GeneralizedBlackScholesProcess> process_;
		Size timeSteps_;
	};


	// template definitions

	template <class T>
	void BinomialVanillaEngine_2<T>::calculate() const {

		DayCounter rfdc = process_->riskFreeRate()->dayCounter();
		DayCounter divdc = process_->dividendYield()->dayCounter();
		DayCounter voldc = process_->blackVolatility()->dayCounter();
		Calendar volcal = process_->blackVolatility()->calendar();

		Real s0 = process_->stateVariable()->value();
		QL_REQUIRE(s0 > 0.0, "negative or null underlying given");
		Volatility v = process_->blackVolatility()->blackVol(
			arguments_.exercise->lastDate(), s0);
		Date maturityDate = arguments_.exercise->lastDate();
		Rate r = process_->riskFreeRate()->zeroRate(maturityDate,
			rfdc, Continuous, NoFrequency);
		Rate q = process_->dividendYield()->zeroRate(maturityDate,
			divdc, Continuous, NoFrequency);
		Date referenceDate = process_->riskFreeRate()->referenceDate();

		// binomial trees with constant coefficient
		Handle<YieldTermStructure> flatRiskFree(
			boost::shared_ptr<YieldTermStructure>(
				new FlatForward(referenceDate, r, rfdc)));
		Handle<YieldTermStructure> flatDividends(
			boost::shared_ptr<YieldTermStructure>(
				new FlatForward(referenceDate, q, divdc)));
		Handle<BlackVolTermStructure> flatVol(
			boost::shared_ptr<BlackVolTermStructure>(
				new BlackConstantVol(referenceDate, volcal, v, voldc)));

		boost::shared_ptr<PlainVanillaPayoff> payoff =
			boost::dynamic_pointer_cast<PlainVanillaPayoff>(arguments_.payoff);
		QL_REQUIRE(payoff, "non-plain payoff given");

		Time maturity = rfdc.yearFraction(referenceDate, maturityDate);

		boost::shared_ptr<StochasticProcess1D> bs(
			new GeneralizedBlackScholesProcess(
				process_->stateVariable(),
				flatDividends, flatRiskFree, flatVol));

		TimeGrid grid(maturity, timeSteps_);

		boost::shared_ptr<T> tree(new T(bs, maturity, timeSteps_,
			payoff->strike()));

		boost::shared_ptr<BlackScholesLattice<T> > lattice(
			new BlackScholesLattice<T>(tree, r, maturity, timeSteps_));

		DiscretizedVanillaOption option(arguments_, *process_, grid);

		option.initialize(lattice, maturity);

		
		/////////////// begin of changes //////////////////
		option.rollback(grid[0]);
		Array value0(option.values());
		
		QL_ENSURE(value0.size() == 3, "Expect 3 nodes in grid at step 0");
		Real p0u = value0[2]; // up
		Real p0 = value0[1];  // mid
		Real p0d = value0[0]; // down
		
		Real s0u = lattice->underlying(0, 2); // up
		Real s0m = lattice->underlying(0, 1); // mid
		Real s0d = lattice->underlying(0, 0); // down
		
		// calculate delta and gamma

		Real delta = ((s0d - s0m) / ((s0u - s0m) * (s0d - s0u)))* p0u + ((s0m - s0u) / ((s0d - s0m) * (s0d - s0u)))* p0d - ((s0u + s0d - 2 * s0m) / ((s0u - s0m) * (s0d - s0m)))* p0;
		Real gamma = 2 * ((s0d - s0u) * p0 - (s0d - s0m) * p0u + (s0u - s0m) * p0d) / ((s0u - s0m)*(s0d - s0u)*(s0d - s0m));

		/////////////// end of changes //////////////////

		// Store results
		results_.value = p0;
		results_.delta = delta;
		results_.gamma = gamma;
		results_.theta = blackScholesTheta(process_,
			results_.value,
			results_.delta,
			results_.gamma);
	}

}


#endif
