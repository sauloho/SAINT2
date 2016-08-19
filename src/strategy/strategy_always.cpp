
#include <iostream>
#include <cassert>
#include "common.h"
#include "strategy_always.h"

// static data members
const char *Strategy_Always::c_type = "always";

Strategy_Always::Strategy_Always()
{
}

Strategy_Always::~Strategy_Always()
{
}

bool Strategy_Always::parse_parameter(const std::string &name,
    const std::string &value)
{
	return false;
}

void Strategy_Always::verify_parameters()
{
}

int Strategy_Always::select(double old_score, const Double_Vec &new_score)
{
	return 0;
}

void Strategy_Always::start_run(Runner *runner)
{
}

void Strategy_Always::end_run(Runner *runner)
{
}

bool Strategy_Always::stop()
{
	// stopping is handled by the Runner
	return false;
}

void Strategy_Always::print_template(std::ostream &out,
	bool commented /*= true*/)
{
	// this is a hidden strategy type
}

