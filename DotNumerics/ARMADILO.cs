
//#define SingleThread

using System;
using System.Collections.Generic;
using System.Text;
using DotNumerics.LinearAlgebra;
using System.Runtime.InteropServices;
using System.Runtime.ExceptionServices;
using System.Security;
using System.Threading.Tasks;
using System.IO;

namespace DotNumerics
{
    public class ARMADILO
    {
        static bool _enabled = true;
        public static bool Enabled
        {
            get
            {
                if (!_enabled) return false;
                return Check();
            }
            set { _enabled = value; }
        }
        public static string dir
        {
            set
            {
                SetDllDirectory(value);
                Enabled = Check();
            }
        }
        static bool? _check = null;
        public static bool Check()
        {
            if (_check == null)
                ReCheck();
            return _check.Value;
        }
        public static bool ReCheck()
        {
            try
            {
                version();
                return (_check = true).Value;
            }
            catch (Exception ex)
            {
                ex.ToString();
                return (_check = false).Value;
            }
        }
        [DllImport(@"armadillo4cs.dll")]
        extern static void fft(int n, double[] v, double[] c_real, double[] c_imag);

        [SecurityCritical, HandleProcessCorruptedStateExceptions]
        public static ComplexVector FFT(Vector b, int N = 0)
        {
            if (N == 0)
                N = b.Count;
            var res = new ComplexVector(N);
            var re = new double[N];
            var im = new double[N];
#if SingleThread
            lock (lock_obj)
#endif
            fft(N, b.data, re, im);
            Parallel.For(0, N, (int i) =>
            {
                res[i] = new Complex(re[i], im[i]);
            });
            return res;
        }
        public static void FFT(Vector b, out Vector result_real, out Vector result_imag, int N = 0)
        {
            if (N == 0)
                N = b.Count;
            var re = new double[N];
            var im = new double[N];
#if SingleThread
            lock (lock_obj)
#endif
            fft(N, b.data, re, im);
            result_real = new Vector(VectorType.Column, false, re);
            result_imag = new Vector(VectorType.Column, false, im);
        }

        /// <summary>
        /// 	equilibrate the system before solving   (matrix A must be square)
        /// </summary>
        public static bool SolverOpt_equilibrate = false;
        /// <summary>
        /// apply iterative refinement to improve solution quality   (matrix A must be square)
        /// </summary>
        public static bool SolverOpt_refine = false;
        /// <summary>
        /// 	fast mode: disable determining solution quality via rcond, disable iterative refinement, disable equilibration
        /// </summary>
        public static bool SolverOpt_fast = false;
        /// <summary>
        /// do not find approximate solutions for rank deficient systems
        /// </summary>
        public static bool SolverOpt_no_approx = true;

        [DllImport("kernel32.dll", CharSet = CharSet.Auto, SetLastError = true)]
        public static extern bool SetDllDirectory(string lpPathName);

        [DllImport(@"armadillo4cs.dll")]
        public extern static double version();

        [DllImport(@"armadillo4cs.dll")]
        extern static void build_time(byte[] str);

        [DllImport(@"armadillo4cs.dll")]
        extern static void build_date(byte[] str);

        [DllImport(@"armadillo4cs.dll")]
        public extern static void set_solve_opt(bool equilibrate, bool refine, bool fast, bool no_approx);
        [DllImport(@"armadillo4cs.dll")]
        public extern static void solve(int n, double[] A, double[] b, double[] x);
        [DllImport(@"armadillo4cs.dll")]
        public extern static void solve_mat(int n, int m, double[] A, double[] b, double[] x);
        [DllImport(@"armadillo4cs.dll")]
        public extern static void cx_solve(int n, double[] A_r, double[] A_i, double[] b_r, double[] b_i, double[] x_r, double[] x_i);
        [DllImport(@"armadillo4cs.dll")]
        public extern static void cx_solve_mat(int n, int m, double[] A_r, double[] A_i, double[] b_r, double[] b_i, double[] x_r, double[] x_i);
        [DllImport(@"armadillo4cs.dll")]
        public extern static void cx_solve2(int n, double[] A_r, double[] A_i, double[] b_r, double[] x_r, double[] x_i);
        /// <summary>
        /// A X = λ B X
        /// </summary> 
        [DllImport(@"armadillo4cs.dll")]
        public extern static void eig_pair(int n, double[] A, double[] B, double[] val_real, double[] val_imag, double[] vec_real, double[] vec_imag);
        /// <summary>
        /// A X = λ B X
        /// </summary> 
        [DllImport(@"armadillo4cs.dll")]
        extern static void cx_eig_pair(int n, double[] A_r, double[] A_i, double[] B_r, double[] B_i, double[] val_real, double[] val_imag, double[] vec_real, double[] vec_imag);
        [DllImport(@"armadillo4cs.dll")]
        public extern static void eig_gen(int n, double[] A, double[] val_real, double[] val_imag, double[] vec_real, double[] vec_imag);

        [DllImport(@"armadillo4cs.dll")]
        public extern static void eigs_gen(int n, double[] A, double[] B, double[] val_real, double[] val_imag, double[] vec_real, double[] vec_imag, int N);

        [DllImport(@"armadillo4cs.dll")]
        public extern static void eig_sym(int n, double[] A, double[] B, double[] val, double[] vec);

        [DllImport(@"armadillo4cs.dll")]
        public extern static void inv(int n, double[] A, double[] Ainv);
        [DllImport(@"armadillo4cs.dll")]
        public extern static void cx_inv(int n, double[] A_r, double[] A_i, double[] Ainv_r, double[] Ainv_i);

        [DllImport(@"armadillo4cs.dll")]
        public extern static void pinv(int n, int m, double[] A, double[] Ainv);

        [DllImport(@"armadillo4cs.dll")]
        public extern static void cx_pinv(int n, int m, double[] A_r, double[] A_i, double[] Ainv_r, double[] Ainv_i);


        [DllImport(@"armadillo4cs.dll")]
        public extern static double det(int n, double[] A);

        [DllImport(@"armadillo4cs.dll")]
        public extern static void cx_det(int n, double[] A_r, double[] A_i, double[] res);


        [DllImport(@"armadillo4cs.dll")]
        public extern static void svd(int n, int m, double[] A, double[] S, double[] U, double[] VT);

        [DllImport(@"armadillo4cs.dll")]
        extern static void mat_multiply(int n, int m, int p, double[] A, double[] B, double[] AB);
        [DllImport(@"armadillo4cs.dll")]
        extern static void cx_mat_multiply_ab(int n, int m, int p, double[] A_r, double[] A_i, double[] B_r, double[] B_i, double[] AB_r, double[] AB_i);
        [DllImport(@"armadillo4cs.dll")]
        extern static void cx_mat_multiply_a(int n, int m, int p, double[] A_r, double[] A_i, double[] B_r, double[] AB_r, double[] AB_i);
        [DllImport(@"armadillo4cs.dll")]
        extern static void cx_mat_multiply_b(int n, int m, int p, double[] A_r, double[] B_r, double[] B_i, double[] AB_r, double[] AB_i);

        [DllImport(@"armadillo4cs.dll")]
        extern static void mv_multiply(int n, int m, double[] A, double[] b, double[] Ab);
        [DllImport(@"armadillo4cs.dll")]
        extern static void mv_multiply1(int n, int m, double[] A_r, double[] A_i, double[] b_r, double[] b_i, double[] Ab_r, double[] Ab_i);
        [DllImport(@"armadillo4cs.dll")]
        extern static void mv_multiply2(int n, int m, double[] A_r, double[] A_i, double[] b, double[] Ab_r, double[] Ab_i);
        [DllImport(@"armadillo4cs.dll")]
        extern static void mv_multiply3(int n, int m, double[] A, double[] b_r, double[] b_i, double[] Ab_r, double[] Ab_i);
        [DllImport(@"armadillo4cs.dll")]
        extern static void roots(int n, double[] P_r, double[] P_i, double[] R_r, double[] R_i);


        private static readonly object lock_obj = new object();


        [SecurityCritical, HandleProcessCorruptedStateExceptions]
        public static DateTime BuildTime()
        {
            var buf = new byte[20];
            build_time(buf);
            var t = Encoding.ASCII.GetString(buf).Replace("\0", "");
            buf = new byte[25];
            build_date(buf);
            var d = Encoding.ASCII.GetString(buf).Replace("\0", "");
            return DateTime.Parse(d + " " + t);
        }
        /// <summary>
        /// A X = λ B X
        /// </summary> 
        [SecurityCritical, HandleProcessCorruptedStateExceptions]
        public static ComplexVector eig_pair(Matrix A, Matrix B, out ComplexMatrix vectors)
        {
            var n = A.RowCount;
            var vr = new double[n];
            var vi = new double[n];
            var vcr = new double[n * n];
            var vci = new double[n * n];
#if SingleThread
            lock (lock_obj)
#endif
            eig_pair(n, A.Data, B.Data, vr, vi, vcr, vci);
            var res = new ComplexVector(n);
            vectors = new ComplexMatrix(n, n);
            for (int i = 0; i < n; i++)
            {
                if (double.IsPositiveInfinity(vr[i]))
                {
                    vr[i] = double.MaxValue;
                    if (double.IsNaN(vi[i]))
                        vi[i] = 0;
                }
                if (double.IsNegativeInfinity(vr[i]))
                {
                    vr[i] = double.MinValue;
                    if (double.IsNaN(vi[i]))
                        vi[i] = 0;
                }
                res[i] = new DotNumerics.Complex(vr[i], vi[i]);
                for (int j = 0; j < n; j++)
                    vectors[i, j] = new DotNumerics.Complex(vcr[i + j * n], vci[i + j * n]);
            }
            return res;
        }
        /// <summary>
        /// A X = λ B X
        /// </summary> 
        [SecurityCritical, HandleProcessCorruptedStateExceptions]
        public static ComplexVector eig_pair(ComplexMatrix A, ComplexMatrix B, out ComplexMatrix vectors)
        {
            var n = A.RowCount;
            var vr = new double[n];
            var vi = new double[n];
            var vcr = new double[n * n];
            var vci = new double[n * n];
#if SingleThread
            lock (lock_obj)
#endif
            cx_eig_pair(n, A.Data_r, A.Data_i, B.Data_r, B.Data_i, vr, vi, vcr, vci);
            var res = new ComplexVector(n);
            vectors = new ComplexMatrix(n, n);
            for (int i = 0; i < n; i++)
            {
                if (double.IsPositiveInfinity(vr[i]))
                {
                    vr[i] = double.MaxValue;
                    if (double.IsNaN(vi[i]))
                        vi[i] = 0;
                }
                if (double.IsNegativeInfinity(vr[i]))
                {
                    vr[i] = double.MinValue;
                    if (double.IsNaN(vi[i]))
                        vi[i] = 0;
                }
                res[i] = new DotNumerics.Complex(vr[i], vi[i]);
                for (int j = 0; j < n; j++)
                    vectors[i, j] = new DotNumerics.Complex(vcr[i + j * n], vci[i + j * n]);
            }
            return res;
        }
        /// <summary>
        /// A X = λ X
        /// </summary> 
        [SecurityCritical, HandleProcessCorruptedStateExceptions]
        public static ComplexVector eig_gen(Matrix A, out ComplexMatrix vectors)
        {
            var n = A.RowCount;
            var vr = new double[n];
            var vi = new double[n];
            var vcr = new double[n * n];
            var vci = new double[n * n];
#if SingleThread
            lock (lock_obj)
#endif
            eig_gen(n, A.Data, vr, vi, vcr, vci);
            var res = new ComplexVector(n);
            vectors = new ComplexMatrix(n, n);
            for (int i = 0; i < n; i++)
            {
                res[i] = new DotNumerics.Complex(vr[i], vi[i]);
                for (int j = 0; j < n; j++)
                    vectors[i, j] = new DotNumerics.Complex(vcr[i + j * n], vci[i + j * n]);
            }
            return res;
        }
        /// <summary>
        /// A X = λ B X
        /// </summary> 
        [SecurityCritical, HandleProcessCorruptedStateExceptions]
        public static ComplexVector eigs_gen(Matrix A, Matrix B, out ComplexMatrix vectors, int N)
        {
            var n = A.RowCount;
            var vr = new double[n];
            var vi = new double[n];
            var vcr = new double[n * n];
            var vci = new double[n * n];
#if SingleThread
            lock (lock_obj)
#endif
            eigs_gen(n, A.Data, B.Data, vr, vi, vcr, vci, N);
            var res = new ComplexVector(n);
            vectors = new ComplexMatrix(n, n);
            for (int i = 0; i < n; i++)
            {
                res[i] = new DotNumerics.Complex(vr[i], vi[i]);
                for (int j = 0; j < n; j++)
                    vectors[i, j] = new DotNumerics.Complex(vcr[i + j * n], vci[i + j * n]);
            }
            return res;
        }
        [SecurityCritical, HandleProcessCorruptedStateExceptions]
        public static Vector eig_sym(Matrix A, Matrix B, out Matrix vectors)
        {
            var n = A.RowCount;
            var res = new Vector(n);
            vectors = new Matrix(n, n);
#if SingleThread
            lock (lock_obj)
#endif
            eig_sym(n, A.Data, B.Data, res.data, vectors.Data);
            return res;
        }
        [SecurityCritical, HandleProcessCorruptedStateExceptions]
        public static Vector solve(Matrix A, Vector b)
        {
            set_solve_opt(SolverOpt_equilibrate, SolverOpt_refine, SolverOpt_fast, SolverOpt_no_approx);
            var x = new Vector(A.RowCount);
            var x_data = x.data;
#if SingleThread
            lock (lock_obj)
#endif
            solve(A.RowCount, A.Data, b.data, x_data);
            return x;
        }
        [SecurityCritical, HandleProcessCorruptedStateExceptions]
        public static Matrix solve(Matrix A, Matrix B)
        {
            set_solve_opt(SolverOpt_equilibrate, SolverOpt_refine, SolverOpt_fast, SolverOpt_no_approx);
            var x = new Vector(A.RowCount, B.ColumnCount);
            var x_data = x.data;
#if SingleThread
            lock (lock_obj)
#endif
            solve_mat(A.RowCount, B.ColumnCount, A.Data, B.Data, x_data);
            return x;
        }
        [SecurityCritical, HandleProcessCorruptedStateExceptions]
        public static Matrix pinv(Matrix A)
        {
            var Ainv = new Matrix(A.ColumnCount, A.RowCount);

#if SingleThread
            lock (lock_obj)
#endif
            pinv(A.RowCount, A.ColumnCount, A.Data, Ainv.Data);
            return Ainv;
        }
        [SecurityCritical, HandleProcessCorruptedStateExceptions]
        public static Matrix inv(Matrix A)
        {
            var Ainv = new Matrix(A.RowCount, A.RowCount);

#if SingleThread
            lock (lock_obj)
#endif
            inv(A.RowCount, A.Data, Ainv.Data);
            return Ainv;
        }

        /// <summary>
        /// P[0] * X^N + ... + P[N],  
        /// P.Count is (N+1), res.Count is N
        /// </summary>
        /// <param name="P"></param>
        /// <returns></returns>
        public static ComplexVector Roots(ComplexVector P)
        {
            var N = P.Count - 1;
            var x_r = new double[N];
            var x_i = new double[N];
            roots(N, P.Data_r, P.Data_i, x_r, x_i);
            var R = new ComplexVector(x_r, x_i);
            return R;
        }

        public static void SetSolveOpt()
        {
            set_solve_opt(SolverOpt_equilibrate, SolverOpt_refine, SolverOpt_fast, SolverOpt_no_approx);
        }

        [SecurityCritical, HandleProcessCorruptedStateExceptions]
        public static ComplexVector solve(ComplexMatrix A, ComplexVector b)
        {
            SetSolveOpt();
            var n = b.Count;
            var x_r = new double[n];
            var x_i = new double[n];
#if SingleThread
            lock (lock_obj)
#endif
            cx_solve(n, A.Data_r, A.Data_i, b.Data_r, b.Data_i, x_r, x_i);
            var x = new ComplexVector(x_r, x_i);
            return x;
        }
        [SecurityCritical, HandleProcessCorruptedStateExceptions]
        public static ComplexVector solve(ComplexMatrix A, Vector b)
        {
            SetSolveOpt();
            var n = b.Count;
            var x_r = new double[n];
            var x_i = new double[n];
#if SingleThread
            lock (lock_obj)
#endif
            cx_solve2(n, A.Data_r, A.Data_i, b.Data, x_r, x_i);
            var x = new ComplexVector(x_r, x_i);
            return x;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A_real"></param>
        /// <param name="A_imag"></param>
        /// <param name="b"></param>
        /// <returns>vector[] {x_real, x_imag}</returns>
        [SecurityCritical, HandleProcessCorruptedStateExceptions]
        public static Vector[] solve(Matrix A_real, Matrix A_imag, Vector b)
        {
            SetSolveOpt();
            var n = b.Count;
            var x_r = new Vector(n);
            var x_i = new Vector(n);
#if SingleThread
            lock (lock_obj)
#endif
            cx_solve2(n, A_real.Data, A_imag.Data, b.Data, x_r.Data, x_i.Data);
            var x = new Vector[] { x_r, x_i };
            return x;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A_real"></param>
        /// <param name="A_imag"></param> 
        /// <returns>vector[] {x_real, x_imag}</returns>
        [SecurityCritical, HandleProcessCorruptedStateExceptions]
        public static Vector[] solve(Matrix A_real, Matrix A_imag, Vector b_real, Vector b_imag)
        {
            SetSolveOpt();
            var n = b_imag.Count;
            var x_r = new Vector(n);
            var x_i = new Vector(n);
#if SingleThread
            lock (lock_obj)
#endif
            cx_solve(n, A_real.Data, A_imag.Data, b_real.Data, b_imag.Data, x_r.Data, x_i.Data);
            var x = new Vector[] { x_r, x_i };
            return x;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A_real"></param>
        /// <param name="A_imag"></param> 
        /// <returns>matrix[] {X_real, X_imag}</returns>
        [SecurityCritical, HandleProcessCorruptedStateExceptions]
        public static Matrix[] solve(Matrix A_real, Matrix A_imag, Matrix b_real, Matrix b_imag)
        {
            SetSolveOpt();
            var n = A_real.RowCount;
            var m = b_real.ColumnCount;
            var x_r = new Matrix(n, m);
            var x_i = new Matrix(n, m);
#if SingleThread
            lock (lock_obj)
#endif
            cx_solve_mat(n, m, A_real.Data, A_imag.Data, b_real.Data, b_imag.Data, x_r.Data, x_i.Data);
            var x = new Matrix[] { x_r, x_i };
            return x;
        }
        [SecurityCritical, HandleProcessCorruptedStateExceptions]
        public static ComplexMatrix inv(ComplexMatrix A)
        {
            var n = A.RowCount;
            var x_r = new double[n * n];
            var x_i = new double[n * n];
#if SingleThread
            lock (lock_obj)
#endif
            cx_inv(A.RowCount, A.Data_r, A.Data_i, x_r, x_i);

            var x = new ComplexMatrix(n, n, x_r, x_i);
            return x;
        }
        [SecurityCritical, HandleProcessCorruptedStateExceptions]
        public static ComplexMatrix pinv(ComplexMatrix A)
        {
            var n = A.RowCount;
            var m = A.ColumnCount;
            var x_r = new double[n * m];
            var x_i = new double[n * m];
#if SingleThread
            lock (lock_obj)
#endif
            cx_pinv(A.RowCount, A.ColumnCount, A.Data_r, A.Data_i, x_r, x_i);

            var x = new ComplexMatrix(m, n, x_r, x_i);
            return x;
        }

        [SecurityCritical, HandleProcessCorruptedStateExceptions]
        public static double det(Matrix A)
        {
#if SingleThread
            lock (lock_obj)
#endif
            return det(A.RowCount, A.Data);
        }
        [SecurityCritical, HandleProcessCorruptedStateExceptions]
        public static Complex det(ComplexMatrix A)
        {
            var res = new double[2];
#if SingleThread
            lock (lock_obj)
#endif
            cx_det(A.RowCount, A.Data_r, A.Data_i, res);
            return new Complex(res[0], res[1]);
        }

        [SecurityCritical, HandleProcessCorruptedStateExceptions]
        public static void svd(Matrix A, out Matrix S, out Matrix U, out Matrix VT)
        {
            int n = A.RowCount;
            int m = A.ColumnCount;
            U = new Matrix(n, n);
            VT = new Matrix(m, m);
            var s_data = new double[Math.Min(m, n)];

#if SingleThread
            lock (lock_obj)
#endif
            svd(n, m, A.Data, s_data, U.Data, VT.Data);
            S = new Matrix(n, m);
            for (int i = 0; i < s_data.Length; i++)
                S[i, i] = s_data[i];
        }
        [SecurityCritical, HandleProcessCorruptedStateExceptions]
        public static Matrix mat_multiply(Matrix A, Matrix B)
        {
            var AB = new Matrix(A.RowCount, B.ColumnCount);

#if SingleThread
            lock (lock_obj)
#endif
            mat_multiply(A.RowCount, A.ColumnCount, B.ColumnCount, A.Data, B.Data, AB.Data);

            return AB;
        }
        [SecurityCritical, HandleProcessCorruptedStateExceptions]
        public static ComplexMatrix mat_multiply(ComplexMatrix A, ComplexMatrix B)
        {
            var n = A.RowCount;
            var m = B.ColumnCount;
            var x_r = new double[n * m];
            var x_i = new double[n * m];
#if SingleThread
            lock (lock_obj)
#endif
            cx_mat_multiply_ab(A.RowCount, A.ColumnCount, B.ColumnCount, A.Data_r, A.Data_i, B.Data_r, B.Data_i, x_r, x_i);

            var x = new ComplexMatrix(n, m, x_r, x_i);
            return x;
        }
        [SecurityCritical, HandleProcessCorruptedStateExceptions]
        public static ComplexMatrix mat_multiply(ComplexMatrix A, Matrix B)
        {
            var n = A.RowCount;
            var m = B.ColumnCount;
            var x_r = new double[n * m];
            var x_i = new double[n * m];
#if SingleThread
            lock (lock_obj)
#endif
            cx_mat_multiply_a(A.RowCount, A.ColumnCount, B.ColumnCount, A.Data_r, A.Data_i, B.Data, x_r, x_i);

            var x = new ComplexMatrix(n, m, x_r, x_i);
            return x;
        }
        [SecurityCritical, HandleProcessCorruptedStateExceptions]
        public static ComplexMatrix mat_multiply(Matrix A, ComplexMatrix B)
        {
            var n = A.RowCount;
            var m = B.ColumnCount;
            var x_r = new double[n * m];
            var x_i = new double[n * m];
#if SingleThread
            lock (lock_obj)
#endif
            cx_mat_multiply_b(A.RowCount, A.ColumnCount, B.ColumnCount, A.Data, B.Data_r, B.Data_i, x_r, x_i);

            var x = new ComplexMatrix(n, m, x_r, x_i);
            return x;
        }
        [SecurityCritical, HandleProcessCorruptedStateExceptions]
        public static Vector mv_multiply(Matrix A, Vector B)
        {
            var x = new Vector(A.RowCount);
            var x_data = x.data;
#if SingleThread
            lock (lock_obj)
#endif
            mv_multiply(A.RowCount, A.ColumnCount, A.Data, B.Data, x_data);

            return x;
        }
        [SecurityCritical, HandleProcessCorruptedStateExceptions]
        public static ComplexVector mv_multiply(ComplexMatrix A, ComplexVector b)
        {
            var n = A.RowCount;
            var m = A.ColumnCount;
            var x_r = new double[n];
            var x_i = new double[n];
#if SingleThread
            lock (lock_obj)
#endif
            mv_multiply1(n, m, A.Data_r, A.Data_i, b.Data_r, b.Data_i, x_r, x_i);

            var x = new ComplexVector(x_r, x_i);
            return x;
        }
        [SecurityCritical, HandleProcessCorruptedStateExceptions]
        public static ComplexVector mv_multiply(ComplexMatrix A, Vector b)
        {
            var n = A.RowCount;
            var m = A.ColumnCount;
            var x_r = new double[n];
            var x_i = new double[n];
#if SingleThread
            lock (lock_obj)
#endif
            mv_multiply2(n, m, A.Data_r, A.Data_i, b.data, x_r, x_i);

            var x = new ComplexVector(x_r, x_i);
            return x;
        }
        [SecurityCritical, HandleProcessCorruptedStateExceptions]
        public static ComplexVector mv_multiply(Matrix A, ComplexVector b)
        {
            var n = A.RowCount;
            var m = A.ColumnCount;
            var x_r = new double[n];
            var x_i = new double[n];
#if SingleThread
            lock (lock_obj)
#endif
            mv_multiply3(n, m, A.Data, b.Data_r, b.Data_i, x_r, x_i);

            var x = new ComplexVector(x_r, x_i);
            return x;
        }


        [DllImport("kernel32.dll", CharSet = CharSet.Unicode, SetLastError = true)]
        static extern int AddDllDirectory(string NewDirectory);
        [DllImport("Kernel32.dll", CharSet = CharSet.Unicode, SetLastError = true)]
        static extern IntPtr LoadLibrary(string path);
        [HandleProcessCorruptedStateExceptions]
        public static void CheckErrorIf(string dir = null)
        {
            dir = dir ?? Directory.GetCurrentDirectory();
            AddDllDirectory(dir);
            LoadLibrary(dir + "\\blas_win64_MT.dll");
            var err = Marshal.GetLastWin32Error() + ",";
            LoadLibrary(dir + "\\lapack_win64_MT.dll");
            err += Marshal.GetLastWin32Error() + ",";
            LoadLibrary(dir + "\\libopenblas.dll");
            err += Marshal.GetLastWin32Error() + ",";
            LoadLibrary(dir + "\\armadillo4cs.dll");
            err += Marshal.GetLastWin32Error() + ",";

            if (!ReCheck())
                throw new Exception("ARMADILO error: " + err);
        }

    }

}