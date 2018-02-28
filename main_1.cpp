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

int main_1() {

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


		boost::shared_ptr<Exercise> europeanExercise(
			new EuropeanExercise(maturity));

		
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
		
		// reference options
		VanillaOption EuropeanOption_blackScholes(payoff, europeanExercise);



		// declare of options 
		VanillaOption europeanOption_Tian_old(payoff, europeanExercise);
		
		VanillaOption europeanOption_Tian_new(payoff, europeanExercise);
		

		/////////////////////// reference options //////////////////
		EuropeanOption_blackScholes.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new AnalyticEuropeanEngine(bsmProcess)));


		/////////////////////// Tian method  ///////////////////////
		cout << "Method =                     " << "Tian" << endl;
		Size timeSteps = 1000;

		// calculate the runtime of the old method
		__int64 BeginTime_oldMethod = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
		europeanOption_Tian_old.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<Tian>(bsmProcess, timeSteps)));
		Real NPV_old = europeanOption_Tian_old.NPV();
		Real delta_old = europeanOption_Tian_old.delta();
		Real gamma_old = europeanOption_Tian_old.gamma();
		__int64 EndTime_oldMethod = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
		

	
		// calculate the runtime of the new method
		__int64 BeginTime_newMethod = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
		europeanOption_Tian_new.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine_2<Tian_2>(bsmProcess, timeSteps)));
		Real NPV_new = europeanOption_Tian_new.NPV();
		Real delta_new = europeanOption_Tian_new.delta();
		Real gamma_new = europeanOption_Tian_new.gamma();
		__int64 EndTime_newMethod = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	
		

		// reference values 
		cout << "Black & Scholes Method" << endl;
		cout << "Price calculated with Black and Scholes method = " << std::scientific << EuropeanOption_blackScholes.NPV() << endl;
		cout << "Delta calculated with Black and Scholes method = " << std::scientific << EuropeanOption_blackScholes.delta() << endl;
		cout << "Gamma calculated with Black and Scholes method = " << std::scientific << EuropeanOption_blackScholes.gamma() << endl;
		cout << endl;

		// statistics for the old method
		cout << "Binomial Tree (Tian) Old method " << endl;
		cout << "Old method runtime (ms) = " << EndTime_oldMethod - BeginTime_oldMethod << endl;
		cout << "Option price calculated with the old method = " << std::scientific << NPV_old << endl;
		cout << "Delta calculated with the old method = " << std::scientific << delta_old << endl;
		cout << "R_BS_delta for the old method = " << std::scientific << (delta_old - EuropeanOption_blackScholes.delta()) / EuropeanOption_blackScholes.delta() << endl;
		cout << "Gamma calculated with the old method = " << std::scientific << gamma_old << endl;
		cout << "R_BS_gamma for the old method = " << std::scientific << (gamma_old - EuropeanOption_blackScholes.gamma()) / EuropeanOption_blackScholes.gamma() << endl;
		cout << endl;

		// statistics for the new method
		cout << "Binomial Tree (Tian) New method " << endl;
		cout << "New method runtime (ms) = " << EndTime_newMethod - BeginTime_newMethod << endl;
		cout << "Option price calculated with the new method = " << std::scientific << NPV_new << endl;
		cout << "Delta calculated with the new method = " << std::scientific << delta_new << endl;
		cout << "R_BS_delta for the new method = " << std::scientific <<(delta_new - EuropeanOption_blackScholes.delta()) / EuropeanOption_blackScholes.delta() << endl;
		cout << "Gamma calculated with the new method = " << std::scientific << gamma_new << endl;
		cout << "R_BS_gamma for the new method = " << std::scientific << (gamma_new - EuropeanOption_blackScholes.gamma()) / EuropeanOption_blackScholes.gamma() << endl;


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