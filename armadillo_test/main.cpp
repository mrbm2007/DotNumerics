#pragma once
#include "stdafx.h"

#include <iostream> 
#include <cstdlib>
#include <ctime> 
#include <iostream>
#include <windows.h>
#include <ctime>
#include <stack>

using namespace std;
using namespace arma;

std::stack<clock_t> tictoc_stack;

void tic() {
	tictoc_stack.push(clock());
}

void toc() {
	std::cout << "Time elapsed: "
		<<1000* ((double)(clock() - tictoc_stack.top())) / CLOCKS_PER_SEC << " ms"
		<< std::endl;
	tictoc_stack.pop();
} 
int main()
{
	mat A = randu<mat>(1000, 1000);
	tic();
	mat c = (A * A).eval();
	toc();

	mat B = pinv(A, 0, "std");        // use default tolerance 
	cout << B[1, 1] ;
}

