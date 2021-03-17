#pragma once
#include "stdafx.h"

#include <iostream> 
#include <cstdlib>
#include <ctime>
#include <omp.h> 
using namespace std;
using namespace arma;

#define ARMA_USE_LAPACK 
#define ARMA_USE_BLAS
#define ARMA_USE_OPENMP 


extern "C" {
	__declspec(dllexport) double version()
	{
		return 10.2010;
	}
	__declspec(dllexport) void build_time(char* str)
	{
		strcpy(str, __TIME__);
	}
	__declspec(dllexport) void build_date(char* str)
	{
		strcpy(str, __DATE__);
	}
	solve_opts::opts solve_opt = solve_opts::opts(solve_opts::no_approx);
	__declspec(dllexport) void set_solve_opt(bool equilibrate, bool refine, bool fast, bool no_approx)
	{
		solve_opt.flags = solve_opts::flag_none;
		if (equilibrate)
			solve_opt.flags += solve_opts::flag_equilibrate;
		if (refine)
			solve_opt.flags += solve_opts::flag_refine;
		if (fast)
			solve_opt.flags += solve_opts::flag_fast;
		if (no_approx)
			solve_opt.flags += solve_opts::flag_no_approx;
	}
	__declspec(dllexport) void solve(int n, double* A_, double* b_, double* x_)
	{
		auto A = mat(A_, n, n, false, false);
		auto b = vec(b_, n, false, false);
		vec x = solve(A, b, solve_opt);
#pragma omp parallel for
		for (register int i = 0;i < n;i++)
			x_[i] = x.mem[i];
	}
	__declspec(dllexport) void solve_mat(int n, int m, double* A_, double* b_, double* x_)
	{
		auto A = mat(A_, n, n, false, false);
		auto b = mat(b_, n, m, false, false);
		mat x = solve(A, b, solve_opt);
#pragma omp parallel for
		for (register int i = 0;i < n * m;i++)
			x_[i] = x.mem[i];
	}

	__declspec(dllexport) void cx_solve(int n, double* A_r, double* A_i, double* b_r, double* b_i, double* x_r, double* x_i)
	{
		auto AA = cx_mat(mat(A_r, n, n, false, false), mat(A_i, n, n, false, false));
		auto bb = cx_vec(vec(b_r, n, false, false), vec(b_i, n, false, false));
		cx_vec xx = solve(AA, bb, solve_opt);
#pragma omp parallel for
		for (register int i = 0;i < n;i++)
		{
			x_r[i] = xx.mem[i].real();
			x_i[i] = xx.mem[i].imag();
		}
	}

	__declspec(dllexport) void cx_solve_mat(int n, int m, double* A_r, double* A_i, double* b_r, double* b_i, double* x_r, double* x_i)
	{
		auto AA = cx_mat(mat(A_r, n, n, false, false), mat(A_i, n, n, false, false));
		auto bb = cx_mat(mat(b_r, n, m, false, false), mat(b_i, n, m, false, false));
		cx_mat xx = solve(AA, bb, solve_opt);
#pragma omp parallel for
		for (register int i = 0;i < n * m;i++)
		{
			x_r[i] = xx.mem[i].real();
			x_i[i] = xx.mem[i].imag();
		}
	}

	__declspec(dllexport) void cx_solve2(int n, double* A_r, double* A_i, double* b_r, double* x_r, double* x_i)
	{
		auto AA = cx_mat(mat(A_r, n, n, false, false), mat(A_i, n, n, false, false));
		auto bb = cx_vec(vec(b_r, n, false, false), vec(n, fill::zeros));
		cx_vec xx = solve(AA, bb, solve_opt);
#pragma omp parallel for
		for (register int i = 0;i < n;i++)
		{
			x_r[i] = xx.mem[i].real();
			x_i[i] = xx.mem[i].imag();
		}
	}

	__declspec(dllexport) void cx_solve_mat2(int n, int m, double* A_r, double* A_i, double* b_r, double* x_r, double* x_i)
	{
		auto AA = cx_mat(mat(A_r, n, n, false, false), mat(A_i, n, n, false, false));
		auto bb = cx_mat(mat(b_r, n, m, false, false), mat(n, m, fill::zeros));
		cx_mat xx = solve(AA, bb, solve_opt);
#pragma omp parallel for
		for (register int i = 0;i < n * m;i++)
		{
			x_r[i] = xx.mem[i].real();
			x_i[i] = xx.mem[i].imag();
		}
	}

	__declspec(dllexport) void eig_pair(int n, double* A_, double* B_,
		double* val_real, double* val_imag, double* vec_real, double* vec_imag)
	{
		auto A = mat(A_, n, n, false, false);
		auto B = mat(B_, n, n, false, false);
		auto val = cx_vec(n);
		auto vec = cx_mat(n, n);
		eig_pair(val, vec, A, B);
#pragma omp parallel for
		for (register int i = 0; i < n; i++)
		{
			val_real[i] = val.mem[i].real();
			val_imag[i] = val.mem[i].imag();
		}
#pragma omp parallel for
		for (register int i = 0; i < n * n; i++)
		{
			vec_real[i] = vec.mem[i].real();
			vec_imag[i] = vec.mem[i].imag();
		}
	}

	__declspec(dllexport) void cx_eig_pair(int n, double* A_r, double* A_i, double* B_r, double* B_i,
		double* val_real, double* val_imag, double* vec_real, double* vec_imag)
	{
		auto A = cx_mat(mat(A_r, n, n, false, false), mat(A_i, n, n, false, false));
		auto B = cx_mat(mat(B_r, n, n, false, false), mat(B_i, n, n, false, false)); 
		auto val = cx_vec(n);
		auto vec = cx_mat(n, n);
		// A*eigvec = B*eigvec*diagmat(eigval)
		eig_pair(val, vec, A, B);
#pragma omp parallel for
		for (register int i = 0; i < n; i++)
		{
			val_real[i] = val.mem[i].real();
			val_imag[i] = val.mem[i].imag();
		}
#pragma omp parallel for
		for (register int i = 0; i < n * n; i++)
		{
			vec_real[i] = vec.mem[i].real();
			vec_imag[i] = vec.mem[i].imag();
		}
	}

	__declspec(dllexport) void eig_gen(int n, double* A_,
		double* val_real, double* val_imag, double* vec_real, double* vec_imag)
	{
		auto A = mat(A_, n, n, false, false);
		auto val = cx_vec(n);
		auto vec = cx_mat(n, n);
		eig_gen(val, vec, A);
#pragma omp parallel for
		for (register int i = 0; i < n; i++)
		{
			val_real[i] = val.mem[i].real();
			val_imag[i] = val.mem[i].imag();
		}
#pragma omp parallel for
		for (register int i = 0; i < n * n; i++) {
			vec_real[i] = vec.mem[i].real();
			vec_imag[i] = vec.mem[i].imag();
		}
	}

	__declspec(dllexport) void eigs_gen(int n, double* A_, double* B_,
		double* val_real, double* val_imag, double* vec_real, double* vec_imag, int N)
	{
		auto A = mat(A_, n, n, false, false);
		auto B = mat(B_, n, n, false, false);
		cx_vec val;
		cx_mat vec;
		auto ba = (inv(B) * A).eval();
		auto bb = sp_mat(n, n);
		for (int i = 0; i < n; i++)
			bb[i] = ba[i];
		eigs_gen(val, vec, bb, N);
#pragma omp parallel for
		for (register int i = 0; i < N; i++)
		{
			val_real[i] = val.mem[i].real();
			val_imag[i] = val.mem[i].imag();
		}
#pragma omp parallel for
		for (register int i = 0; i < n * N; i++)
		{
			vec_real[i] = vec.mem[i].real();
			vec_imag[i] = vec.mem[i].imag();
		}
	}

	__declspec(dllexport) void eig_sym(int n, double* A_, double* B_,
		double* val_, double* vec_)
	{
		auto A = mat(A_, n, n, false, false);
		auto B = mat(B_, n, n, false, false);
		auto val = vec(val_, n, false, false);
		auto vec = mat(vec_, n, n, false, false);
		eig_sym(val, vec, inv(B) * A);
#pragma omp parallel for
		for (register int i = 0; i < n; i++)
			val_[i] = val.mem[i];
#pragma omp parallel for
		for (register int i = 0; i < n * n; i++)
			vec_[i] = vec.mem[i];
	}

	__declspec(dllexport) double det(int n, double* A_)
	{
		auto A = mat(A_, n, n, false, false);
		return arma::det(A);
	}
	__declspec(dllexport) void cx_det(int n, double* A_r, double* A_i, double* res)
	{
		auto A = cx_mat(mat(A_r, n, n, false, false), mat(A_i, n, n, false, false));
		auto res_ = arma::det(A);
		res[0] = res_.real();
		res[1] = res_.imag();
	}

	__declspec(dllexport) void roots(int n, double* P_r, double* P_i, double* R_r, double* R_i)
	{
		auto P = cx_vec(vec(P_r, n+1,  false), vec(P_i, n+1, false));
		auto R = roots(P).eval();
#pragma omp parallel for
		for (register int i = 0;i < n;i++) {
			R_r[i] = R.mem[i].real();
			R_i[i] = R.mem[i].imag();
		}
	}

	__declspec(dllexport) void inv(int n, double* A_, double* Ainv_)
	{
		auto A = mat(A_, n, n, false, false);
		auto Ainv = inv(A).eval();
#pragma omp parallel for
		for (register int i = 0;i < n * n;i++)
			Ainv_[i] = Ainv.mem[i];
	}
	__declspec(dllexport) void cx_inv(int n, double* A_r, double* A_i, double* Ainv_r, double* Ainv_i)
	{
		auto A = cx_mat(mat(A_r, n, n, false, false), mat(A_i, n, n, false, false));
		cx_mat Ainv = inv(A).eval();
#pragma omp parallel for
		for (register int i = 0;i < n * n;i++)
		{
			Ainv_r[i] = Ainv.mem[i].real();
			Ainv_i[i] = Ainv.mem[i].imag();
		}
	}
	__declspec(dllexport) void pinv(int n, int m, double* A_, double* Ainv_)
	{
		auto A = mat(A_, n, m, false, false);
		auto Ainv = pinv(A, 0, "std").eval();
#pragma omp parallel for
		for (register int i = 0;i < n * m;i++)
			Ainv_[i] = Ainv.mem[i];
	}
	__declspec(dllexport) void cx_pinv(int n, int m, double* A_r, double* A_i, double* Ainv_r, double* Ainv_i)
	{
		auto A = cx_mat(mat(A_r, n, n, false, false), mat(A_i, n, n, false, false));
		cx_mat Ainv = pinv(A, 0, "std").eval();
#pragma omp parallel for
		for (register int i = 0;i < n * n;i++)
		{
			Ainv_r[i] = Ainv.mem[i].real();
			Ainv_i[i] = Ainv.mem[i].imag();
		}
	}

	__declspec(dllexport) void svd(int n, int m, double* A_, double* s_, double* U_, double* VT_)
	{
		auto A = mat(A_, n, m, false, false);
		auto U = mat(U_, n, n, false, false);
		auto mn = min(n, m);
		auto s = vec(s_, mn, false, false);
		auto V = mat(VT_, m, m, false, false);
		svd(U, s, V, A, "std");
		V = V.t().eval();
#pragma omp parallel for
		for (register int i = 0; i < m * m; i++)
			VT_[i] = V.mem[i];
#pragma omp parallel for
		for (register int i = 0; i < mn * mn; i++)
			s_[i] = s.mem[i];
#pragma omp parallel for
		for (register int i = 0; i < n * n; i++)
			U_[i] = U.mem[i];
	}
	/*
	__declspec(dllexport) void svd(int n, int m,
		double* A_r_, double* A_i_, double* s_r_, double* s_i_,
		double* U_r_, double* U_i_, double* VT_r_, double* VT_i_)
	{
		auto A = cx_mat(mat(A_r_, n, m, false, false), mat(A_i_, n, m, false, false));
		auto U = cx_mat(mat(U_r_, n, n, false, false), mat(U_i_, n, n, false, false));
		auto s = cx_mat(vec(s_r_, min(n, m), false, false), vec(s_i_, min(n, m), false, false));
		auto V = cx_mat(mat(VT_r_, m, m, false, false), mat(VT_i_, m, m, false, false));
		svd(U, s, V, A, "std");
		V = V.t().eval();
	}*/

	__declspec(dllexport) void fft(int v_size, double* v, double* c_real, double* c_imag)
	{
		auto b = vec(v, v_size, false, false);
		cx_vec Y = fft(b, v_size);
#pragma omp parallel for
		for (register int i = 0; i < v_size; i++)
		{
			c_real[i] = Y.mem[i].real();
			c_imag[i] = Y.mem[i].imag();
		}
	}
	__declspec(dllexport) void mat_multiply(int n, int m, int p, double* A_, double* B_, double* AB_)
	{
		/*
#pragma omp parallel for
		for (register int i = 0; i < n; i++)
		{
			auto im = i*m;
			auto ip = i*p;
#pragma omp parallel for
			for (int j = 0; j < p; j++)
			{
				auto register sum = 0.0;
				for (int k = 0; k < m; k++)
					sum += A_[im + k] * B_[k*p + j];
				AB_[ip + j] = sum;
			}
		}
		return;*/
		auto A = mat(A_, n, m, false, false);
		auto B = mat(B_, m, p, false, false);
		auto AB = (A * B).eval();
#pragma omp parallel for
		for (register int i = 0; i < n * p; i++)
			AB_[i] = AB.mem[i];
	}
	__declspec(dllexport) void cx_mat_multiply_ab(int n, int m, int p, double* A_r, double* A_i, double* B_r, double* B_i, double* AB_r, double* AB_i)
	{
		auto A = cx_mat(mat(A_r, n, m, false, false), mat(A_i, n, m, false, false));
		auto B = cx_mat(mat(B_r, m, p, false, false), mat(B_i, m, p, false, false));
		auto AB = (A * B).eval();
#pragma omp parallel for
		for (register int i = 0; i < n * p; i++)
		{
			AB_r[i] = AB.mem[i].real();
			AB_i[i] = AB.mem[i].imag();
		}
	}
	__declspec(dllexport) void cx_mat_multiply_a(int n, int m, int p, double* A_r, double* A_i, double* B_r, double* AB_r, double* AB_i)
	{
		auto A = cx_mat(mat(A_r, n, m, false, false), mat(A_i, n, m, false, false));
		auto B = mat(B_r, m, p, false, false);
		auto AB = (A * B).eval();
#pragma omp parallel for
		for (register int i = 0; i < n * p; i++)
		{
			AB_r[i] = AB.mem[i].real();
			AB_i[i] = AB.mem[i].imag();
		}
	}
	__declspec(dllexport) void cx_mat_multiply_b(int n, int m, int p, double* A_r, double* B_r, double* B_i, double* AB_r, double* AB_i)
	{
		auto A = mat(A_r, n, m, false, false);
		auto B = cx_mat(mat(B_r, m, p, false, false), mat(B_i, m, p, false, false));
		auto AB = (A * B).eval();
#pragma omp parallel for
		for (register int i = 0; i < n * p; i++)
		{
			AB_r[i] = AB.mem[i].real();
			AB_i[i] = AB.mem[i].imag();
		}
	}

	__declspec(dllexport) void mv_multiply(int n, int m, double* A_, double* b_, double* Ab_)
	{
		auto A = mat(A_, n, m, false, false);
		auto b = vec(b_, m, false, false);
		auto Ab = (A * b).eval();
#pragma omp parallel for
		for (register int i = 0; i < n; i++)
			Ab_[i] = Ab.mem[i];
	}
	__declspec(dllexport) void mv_multiply1(int n, int m, double* A_r, double* A_i, double* b_r, double* b_i, double* Ab_r, double* Ab_i)
	{
		auto A = cx_mat(mat(A_r, n, m, false, false), mat(A_i, n, m, false, false));
		auto b = cx_vec(vec(b_r, m, false, false), vec(b_i, m, false, false));
		auto Ab = (A * b).eval();
#pragma omp parallel for
		for (register int i = 0; i < n; i++)
		{
			Ab_r[i] = Ab.mem[i].real();
			Ab_i[i] = Ab.mem[i].imag();
		}
	}
	__declspec(dllexport) void mv_multiply2(int n, int m, double* A_r, double* A_i, double* b_, double* Ab_r, double* Ab_i)
	{
		auto A = cx_mat(mat(A_r, n, m, false, false), mat(A_i, n, m, false, false));
		auto b = vec(b_, m, false, false);
		auto Ab = (A * b).eval();
#pragma omp parallel for
		for (register int i = 0; i < n; i++)
		{
			Ab_r[i] = Ab.mem[i].real();
			Ab_i[i] = Ab.mem[i].imag();
		}
	}
	__declspec(dllexport) void mv_multiply3(int n, int m, double* A_, double* b_r, double* b_i, double* Ab_r, double* Ab_i)
	{
		auto A = mat(A_, n, m, false, false);
		auto b = cx_vec(vec(b_r, m, false, false), vec(b_i, m, false, false));
		auto Ab = (A * b).eval();
#pragma omp parallel for
		for (register int i = 0; i < n; i++)
		{
			Ab_r[i] = Ab.mem[i].real();
			Ab_i[i] = Ab.mem[i].imag();
		}
	}
}
