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

int main() {

	try {
		Size widths[] = { 35, 14, 14, 14 };
		boost::timer timer;
		std::cout << std::endl;

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
		std::cout << "Option type =                " << type << std::endl;
		std::cout << "Maturity =                   " << maturity << std::endl;
		std::cout << "Underlying price =           " << underlying << std::endl;
		std::cout << "Strike =                     " << strike << std::endl;
		std::cout << "Risk-free interest rate =    " << io::rate(riskFreeRate) << std::endl;
		std::cout << "Dividend yield =             " << io::rate(dividendYield) << std::endl;
		std::cout << "Volatility =                 " << io::volatility(volatility) << std::endl;

		std::cout << std::endl;
		std::string method;
		std::cout << std::endl;

		//// daclare the 3 types of option-exercice: european, bermudan and american
		//std::vector<Date> exerciseDates;
		//for (Integer i = 1; i <= 4; i++)
		//	exerciseDates.push_back(settlementDate + 3 * i*Months);

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



		// declare 3 types of options for each calculation method
		VanillaOption europeanOption_Tian_old(payoff, europeanExercise);
		
		VanillaOption europeanOption_Tian_new(payoff, europeanExercise);
		

		/////////////////////// reference options //////////////////
		EuropeanOption_blackScholes.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new AnalyticEuropeanEngine(bsmProcess)));


		/////////////////////// Tian method  ///////////////////////
		method = "Binomial Tian";
		Size timeSteps = 1000;

		__int64 BeginTime_oldMethod = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
		europeanOption_Tian_old.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine<Tian>(bsmProcess, timeSteps)));
		__int64 EndTime_oldMethod = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
		

		std::cout << "time taken by the old method   " << EndTime_oldMethod - BeginTime_oldMethod << std::endl;


		__int64 BeginTime_newMethod = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

		europeanOption_Tian_new.setPricingEngine(boost::shared_ptr<PricingEngine>(
			new BinomialVanillaEngine_2<Tian_2>(bsmProcess, timeSteps)));
		__int64 EndTime_newMethod = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	
		std::cout << "time taken by the new method   " << EndTime_newMethod - BeginTime_newMethod << std::endl;


		std::cout << std::setw(widths[0]) << std::left << "Method"
			<< std::setw(widths[1]) << std::left << "European" 
			<< std::endl;
		std::cout << std::setw(widths[0]) << std::left << method
			<< std::fixed
			<< std::setw(widths[1]) << std::left << europeanOption_Tian_old.NPV()
			<< std::endl;

		std::cout << std::setw(widths[0]) << std::left << method
			<< std::fixed
			<< std::setw(widths[1]) << std::left << europeanOption_Tian_new.NPV()
			<< std::endl;

		std::cout << "Delta calculated with the old method = " << europeanOption_Tian_old.delta() << std::endl;
		std::cout << "Delta calculated with the new method = " << europeanOption_Tian_new.delta() << std::endl;
		std::cout << "Delta calculated with Black and Scholes method = " << EuropeanOption_blackScholes.delta() << std::endl;
		std::cout << "The ratio is    " << (europeanOption_Tian_new.delta() - EuropeanOption_blackScholes.delta()) / EuropeanOption_blackScholes.delta() << endl;




		std::cout << "gamma calculated with the old method = " << europeanOption_Tian_old.gamma() << std::endl;
		std::cout << "gamma calculated with the new method = " << europeanOption_Tian_new.gamma() << std::endl;
		std::cout << "gamma calculated with Black and Scholes method = " << EuropeanOption_blackScholes.gamma() << std::endl;
		cout << "The ratio is   " << (europeanOption_Tian_new.gamma() - EuropeanOption_blackScholes.gamma()) / EuropeanOption_blackScholes.gamma() << endl;


		std::cin.ignore();

		return 0;

	}
	catch (std::exception& e) {
		std::cerr << e.what() << std::endl;
		std::cin.ignore();
		return 1;
	}
	catch (...) {
		std::cerr << "unknown error" << std::endl;
		std::cin.ignore();
		return 1;
	}
}