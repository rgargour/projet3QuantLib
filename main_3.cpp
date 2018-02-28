/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
#pragma warning(disable : 4996)

#include <ql/qldefines.hpp>
#ifdef BOOST_MSVC
#  include <ql/auto_link.hpp>
#endif
#include <ql/instruments/vanillaoption.hpp>
#include <ql/pricingengines/vanilla/binomialengine.hpp>
#include <ql/pricingengines/vanilla/fdeuropeanengine.hpp>
#include <ql/pricingengines/vanilla/fdbermudanengine.hpp>
#include <ql/pricingengines/vanilla/fdamericanengine.hpp>
#include <ql/pricingengines/vanilla/mceuropeanengine.hpp>
#include <ql/pricingengines/vanilla/mcamericanengine.hpp>
#include <ql/pricingengines/vanilla/analytichestonengine.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/math/distributions/binomialdistribution.hpp>
#include <ql/stochasticprocess.hpp>
#include <ql/utilities/dataformatters.hpp>
#include <ql/pricingengines/vanilla/baroneadesiwhaleyengine.hpp>
#include <ql/pricingengines/vanilla/bjerksundstenslandengine.hpp>
#include "binomialtree.hpp"
#include "binomialengine.hpp"
#include <ql/methods/lattices/binomialtree.hpp>
#include <ql/pricingengines/vanilla/binomialengine.hpp>
#include <iostream>
#include <boost/timer.hpp>
#include <iomanip>
#include <ql/time/calendars/target.hpp>
#include <ctime> 
#include <chrono>

using std::cout;
using std::endl;

using namespace QuantLib;

int main_3() {

	try {
		Size widths[] = { 35, 14, 14, 14 };
		boost::timer timer;
		cout << endl;

		// set up dates
		Calendar calendar = TARGET();
		Date todaysDate(13, February, 2018);
		Date settlementDate(17, February, 2018);
		Settings::instance().evaluationDate() = todaysDate;

		// our options
		Option::Type type(Option::Call);
		Real underlying = 36;
		Real strike = 40;
		Spread dividendYield = 0.00;
		Rate riskFreeRate = 0.06;
		Volatility volatility = 0.20;
		Date maturity(13, February, 2019);
		DayCounter dayCounter = Actual365Fixed();

		// to display the caracteristics of the option
		cout << "Option type =                " << type << endl;
		cout << "Maturity =                   " << maturity << endl;
		cout << "Underlying price =           " << underlying << endl;
		cout << "Strike =                     " << strike << endl;
		cout << "Risk-free interest rate =    " << io::rate(riskFreeRate) << endl;
		cout << "Dividend yield =             " << io::rate(dividendYield) << endl;
		cout << "Volatility =                 " << io::volatility(volatility) << endl;

		cout << endl;
		std::string method;
		cout << endl;



		boost::shared_ptr<Exercise> americanExercise(
			new AmericanExercise(settlementDate,
				maturity));

		Handle<Quote> underlyingH(
			boost::shared_ptr<Quote>(new SimpleQuote(underlying)));

		// bootstrap the yield/dividend/vol curves
		Handle<YieldTermStructure> flatTermStructure(
			boost::shared_ptr<YieldTermStructure>(
				new FlatForward(settlementDate, riskFreeRate, dayCounter)));
		Handle<YieldTermStructure> flatDividendTS(
			boost::shared_ptr<YieldTermStructure>(
				new FlatForward(settlementDate, dividendYield, dayCounter)));
		Handle<BlackVolTermStructure> flatVolTS(
			boost::shared_ptr<BlackVolTermStructure>(
				new BlackConstantVol(settlementDate, calendar, volatility,
					dayCounter)));
		boost::shared_ptr<StrikedTypePayoff> payoff(
			new PlainVanillaPayoff(type, strike));
		boost::shared_ptr<BlackScholesMertonProcess> bsmProcess(
			new BlackScholesMertonProcess(underlyingH, flatDividendTS,
				flatTermStructure, flatVolTS));


		// declare options 
		Size timeSteps = 1000;
		// binomial tree tian 
		VanillaOption americanOption_old_Tian(payoff, americanExercise);
		americanOption_old_Tian.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<Tian>(bsmProcess, timeSteps)));

		VanillaOption americanOption_new_Tian(payoff, americanExercise);
		americanOption_new_Tian.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine_2<Tian_2>(bsmProcess, timeSteps)));

		VanillaOption americanOption_BaroneAdesiWhaley(payoff, americanExercise);
		americanOption_BaroneAdesiWhaley.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BaroneAdesiWhaleyApproximationEngine(bsmProcess)));

		VanillaOption americanOption_BjerksundStensland(payoff, americanExercise);
		americanOption_BjerksundStensland.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BjerksundStenslandApproximationEngine(bsmProcess)));

		VanillaOption americanOption_FDA(payoff, americanExercise);
		americanOption_FDA.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new FDAmericanEngine<CrankNicolson>(bsmProcess,
				timeSteps, timeSteps - 1)));

		// statistics for the old method
		cout << "Binomial Tree (Tian) Old method " << endl;
		cout << "Delta  = " << std::scientific << americanOption_old_Tian.delta() << endl;
		cout << "Binomial Tree (Tian) New method " << endl;
		cout << "Delta  = " << std::scientific << americanOption_new_Tian.delta() << endl;
		cout << "Barone-Adesi/Whaley" << endl;
		cout << "Delta  = " << std::scientific << americanOption_BaroneAdesiWhaley.delta() << endl;
		cout << "Bjerksund/Stensland" << endl;
		cout << "Delta  = " << std::scientific << americanOption_BjerksundStensland.delta() << endl;
		cout << "Finite differences" << endl;
		cout << "Delta  = " << std::scientific << americanOption_FDA.delta() << endl;
		cout << endl;

		

		std::cin.ignore();

		return 0;

	}
	catch (std::exception& e) {
		std::cerr << e.what() << endl;
		std::cin.ignore();
		return 1;
	}
	catch (...) {
		std::cerr << "unknown error" << endl;
		std::cin.ignore();
		return 1;
	}
}