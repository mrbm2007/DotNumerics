#region Copyright © 2009, De Santiago-Castillo JA. All rights reserved.

//Copyright © 2009 Jose Antonio De Santiago-Castillo 
//E-mail:JAntonioDeSantiago@gmail.com
//Web: www.DotNumerics.com
//
#endregion

using System;
using System.Collections.Generic;
using System.Text;
using System.ComponentModel;
using System.Drawing;
using DotNumerics.LinearAlgebra.CSLapack;

namespace DotNumerics.LinearAlgebra
{

    /// <summary>
    /// Computes the solution to a system of linear equations.
    /// </summary>
    public sealed class LinearEquations
    {

        #region Fields

        DGESV _dgesv;
        DGBSV _dgbsv;
        DGTSV _dgtsv;

        #endregion


        /// <summary>
        /// Initializes a new instance of the LinearEquations class.
        /// </summary>
        public LinearEquations()
        {

        }

        #region Public LU Solver


        #region General Matrix
        /// <summary>
        /// Computes the solution to a real system of linear equations A * X = B, where A is a general matrix. 
        /// </summary>
        /// <param name="A">The square matrix.</param>
        /// <param name="B">The vector containing the right-hand side of the linear system.</param>
        /// <returns>A vector containing the solution to the linear system of equations.</returns>
        public static Vector Solve(Matrix A, Vector B)
        {
            return new LinearEquations().Solve_(A, B);
        }
        /// <summary>
        /// Computes the solution to a real system of linear equations A * X = B, where A is a general matrix. 
        /// </summary>
        /// <param name="A">The square matrix.</param>
        /// <param name="b">The vector containing the right-hand side of the linear system.</param>
        /// <returns>A vector containing the solution to the linear system of equations.</returns>
        public Vector Solve_(Matrix A, Vector b)
        {
            CheckDimensions(A, b);
            if (ARMADILO.Enabled)
                return ARMADILO.solve(A, b);
            else
            {
                Matrix Solution = b.Clone();
                Matrix AClon = A.Clone();

                this.SolveInplace(AClon, Solution);

                return Solution.GetColumnVector(0);
            }
        }
        /// <summary>
        /// Computes the solution to a real system of linear equations A * X = B, where A is a general matrix. 
        /// for each column in B
        /// </summary>
        /// <param name="A">The square matrix.</param>
        /// <param name="B">The vector containing the right-hand side of the linear system.</param>
        /// <returns>A vector containing the solution to the linear system of equations.</returns>
        public Matrix Solve_(Matrix A, Matrix B)
        {
            CheckDimensions(A, B);
            if (ARMADILO.Enabled)
                return ARMADILO.solve(A, B);
            else
            {
                Matrix Solution = B.Clone();
                Matrix AClon = A.Clone();

                this.SolveInplace(AClon, Solution);

                return Solution;
            }
        }
        /// <summary>
        /// Computes the solution to a real system of linear equations A * X = B, where A is a general matrix. 
        /// </summary>
        /// <param name="A">The square matrix.</param>
        /// <param name="B">The vector containing the right-hand side of the linear system.</param>
        /// <returns>A vector containing the solution to the linear system of equations.</returns>
        public static ComplexVector Solve(ComplexMatrix A, ComplexVector B)
        {
            CheckDimensions(A, B);
            if (ARMADILO.Enabled)
                return ARMADILO.solve(A, B);
            else
                throw new NotImplementedException("ComplexMatrix.Solve: only with armadillo");
        }
        /// <summary>
        /// Computes the solution to a real system of linear equations A * X = B, where A is a general matrix. 
        /// </summary>
        /// <param name="A">The square matrix.</param>
        /// <param name="B">The vector containing the right-hand side of the linear system.</param>
        /// <returns>A vector containing the solution to the linear system of equations.</returns>
        public static ComplexVector Solve(ComplexMatrix A, Vector B)
        {
            CheckDimensions(A, B);
            if (ARMADILO.Enabled)
                return ARMADILO.solve(A, B);
            else
                throw new NotImplementedException("ComplexMatrix.Solve: only with armadillo");
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A_real"></param>
        /// <param name="A_imag"></param>
        /// <param name="B"></param>
        /// <returns>vector[] {x_real, x_imag}</returns>
        public static Vector[] Solve(Matrix A_real, Matrix A_imag, Vector B)
        {
            CheckDimensions(A_real, B);
            CheckDimensions(A_imag, B);
            if (ARMADILO.Enabled)
                return ARMADILO.solve(A_real, A_imag, B);
            else
                throw new NotImplementedException("ComplexMatrix.Solve: only with armadillo");
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A_real"></param>
        /// <param name="A_imag"></param> 
        /// <returns>vector[] {x_real, x_imag}</returns>
        public static Vector[] Solve(Matrix A_real, Matrix A_imag, Vector b_real, Vector b_imag)
        {
            CheckDimensions(A_real, b_real);
            CheckDimensions(A_imag, b_imag);
            if (b_real.Count != b_imag.Count)
                CheckDimensions(A_imag, b_real);// will fail
            if (ARMADILO.Enabled)
                return ARMADILO.solve(A_real, A_imag, b_real, b_imag);
            else
                throw new NotImplementedException("ComplexMatrix.Solve: only with armadillo");
        }
        /// <summary>
        /// solve for each column in b
        /// </summary>
        /// <param name="A_real"></param>
        /// <param name="A_imag"></param> 
        /// <returns>Matrix[] {X_real, X_imag}</returns>
        public static Matrix[] Solve(Matrix A_real, Matrix A_imag, Matrix b_real, Matrix b_imag)
        {
            CheckDimensions(A_real, b_real);
            CheckDimensions(A_imag, b_imag);
            if (b_real.RowCount != b_imag.RowCount)
                CheckDimensions(A_imag, b_real);// will fail
            if (ARMADILO.Enabled)
                return ARMADILO.solve(A_real, A_imag, b_real, b_imag);
            else
                throw new NotImplementedException("ComplexMatrix.Solve: only with armadillo");
        }

        /// <summary>
        /// Computes the solution to a real system of linear equations A * X = B, where A is a general matrix. 
        /// </summary>
        /// <param name="A">The square matrix.</param>
        /// <param name="B">The vector containing the right-hand side of the linear system.</param>
        /// <returns>A vector containing the solution to the linear system of equations.</returns>
        public static double[] Solve(double[,] A, double[] B)
        {

            var Solution = new Vector(B);
            var AClon = new Matrix(A);

            return Solve(AClon, Solution).data;

        }


        /// <summary>
        /// Computes the solution to a real system of linear equations A * X = B, where A is a general matrix. 
        /// </summary>
        /// <param name="A">The square matrix.</param>
        /// <param name="B">The matrix containing the right-hand side of the linear system.</param>
        /// <returns>A matrix containing the solution to the linear system of equations.</returns>
        public Matrix Solve(Matrix A, Matrix B)
        {
            CheckDimensions(A, B);

            Matrix Solution = B.Clone();
            Matrix AClon = A.Clone();

            this.SolveInplace(AClon, Solution);


            return Solution;
        }


        /// <summary>
        /// In place, Computes the solution to a real system of linear equations A * X = B
        /// </summary>
        /// <param name="A">The square matrix.</param>
        /// <param name="B">The matrix containing the right-hand side of the linear system.</param>
        /// <returns>A matrix containing the solution to the linear system of equations.</returns>
        private void SolveInplace(Matrix A, Matrix B)
        {

            if (this._dgesv == null) this._dgesv = new DGESV();

            CheckDimensions(A, B);

            int numberEquations = A.RowCount;
            int numberSolutions = B.ColumnCount;
            double[] AMatrix = A.Data;
            double[] BMatrix = B.Data;

            int[] IPIV = new int[numberEquations];
            int Info = 0;

            this._dgesv.Run(numberEquations, numberSolutions, ref AMatrix, 0, numberEquations, ref IPIV, 0, ref BMatrix, 0, numberEquations, ref Info);


            #region Error
            /// = 0:  successful exit
            /// .LT. 0:  if INFO = -i, the i-th argument had an illegal value
            /// .GT. 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
            /// has been completed, but the factor U is exactly
            /// singular, so the solution could not be computed.

            if (Info < 0)
            {
                string infoSTg = Math.Abs(Info).ToString();
                throw new ArgumentException("the " + infoSTg + " -th argument had an illegal value");
            }
            else if (Info > 0)
            {
                string infoSTg = Math.Abs(Info).ToString();
                throw new Exception("U(" + infoSTg + "," + infoSTg + ") is exactly zero.  The factorization has been completed, but the factor U is exactly singular, so the solution could not be computed.");
            }

            #endregion

        }






        #endregion


        #region BandMatrix


        /// <summary>
        /// Computes the solution to a real system of linear equations
        /// A * X = B, where A is a band matrix.
        /// </summary>
        /// <param name="A">The band matrix.</param>
        /// <param name="B">The vector containing the right-hand side of the linear system.</param>
        /// <returns>A vector containing the solution to the linear system of equations.</returns>
        public Vector Solve(BandMatrix A, Vector B)
        {
            Matrix myB = B;
            Vector solution = this.Solve(A, myB).GetColumnVector(0);
            return this.Solve(A, B);
        }


        /// <summary>
        /// Computes the solution to a real system of linear equations
        /// A * X = B, where A is a band matrix.
        /// </summary>
        /// <param name="A">The band matrix.</param>
        /// <param name="B">The matrix containing the right-hand side of the linear system.</param>
        /// <returns>A matrix containing the solution to the linear system of equations.</returns>
        public Matrix Solve(BandMatrix A, Matrix B)
        {
            if (this._dgbsv == null) this._dgbsv = new DGBSV();

            CheckDimensions(A, B);

            Matrix Solution = B.Clone();
            Matrix ExtendedMatrix = A.GetBandPackedMatrix();
            double[] BAData = ExtendedMatrix.Data;
            double[] SolutionData = Solution.Data;

            int[] IPIV = new int[A.RowCount];
            int Info = 0;

            this._dgbsv.Run(A.RowCount, A.LowerBandWidth, A.UpperBandWidth, B.ColumnCount, ref BAData, 0, ExtendedMatrix.RowCount, ref IPIV, 0, ref SolutionData, 0, Solution.RowCount, ref Info);


            #region Error
            /// = 0:  successful exit
            /// .LT. 0:  if INFO = -i, the i-th argument had an illegal value
            /// .GT. 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
            /// has been completed, but the factor U is exactly
            /// singular, and the solution has not been computed.
            ///</param>

            if (Info < 0)
            {
                string infoSTg = Math.Abs(Info).ToString();
                throw new ArgumentException("the " + infoSTg + " -th argument had an illegal value");
            }
            else if (Info > 0)
            {
                string infoSTg = Math.Abs(Info).ToString();
                throw new Exception("U(" + infoSTg + "," + infoSTg + ") is exactly zero.  The factorization has been completed, but the factor U is exactly singular, so the solution could not be computed.");
            }

            #endregion




            return Solution;
        }


        #endregion


        /// <summary>
        /// Computes the solution to a real system of linear equations
        /// A * X = B, where A is a tridiagonal matrix.
        /// </summary>
        /// <param name="A">The tridiagonal matrix.</param>
        /// <param name="B">The matrix containing the right-hand side of the linear system.</param>
        /// <returns>A matrix containing the solution to the linear system of equations.</returns>
        public Matrix Solve(TridiagonalMatrix A, Matrix B)
        {

            if (this._dgtsv == null) this._dgtsv = new DGTSV();

            CheckDimensions(A, B);

            Matrix Solution = B.Clone();
            double[] SolutionData = Solution.Data;
            double[] Diagonal;
            double[] SubDiagonal;
            double[] SuperDiagonal;
            A.GetPackedMatrix(out SubDiagonal, out SuperDiagonal, out Diagonal);

            int[] IPIV = new int[A.RowCount];
            int Info = 0;

            this._dgtsv.Run(A.RowCount, B.ColumnCount, ref SubDiagonal, 0, ref Diagonal, 0, ref SuperDiagonal, 0, ref SolutionData, 0, Solution.RowCount, ref Info);


            #region Error
            /// = 0: successful exit
            /// .LT. 0: if INFO = -i, the i-th argument had an illegal value
            /// .GT. 0: if INFO = i, U(i,i) is exactly zero, and the solution
            /// has not been computed.  The factorization has not been
            /// completed unless i = N.
            if (Info < 0)
            {
                string infoSTg = Math.Abs(Info).ToString();
                throw new ArgumentException("the " + infoSTg + " -th argument had an illegal value");
            }
            else if (Info > 0)
            {
                string infoSTg = Math.Abs(Info).ToString();
                throw new Exception("U(" + infoSTg + "," + infoSTg + ") is exactly zero.  and the solution has not been computed.  The factorization has not been completed unless i = N.");
            }

            #endregion

            return Solution;
        }



        private static void CheckDimensions(BaseMatrix matrixA, BaseMatrix matrixB)
        {
            if (matrixA.IsSquare != true)
            {
                throw new System.ArgumentException("Matrix A is not a square matrix.");
            }


            if (matrixA.RowCount != matrixB.RowCount)
            {
                throw new System.ArgumentException("Matrix dimensions are not valid.");
            }
        }
        private static void CheckDimensions(ComplexMatrix matrixA, Vector matrixB)
        {
            if (matrixA.IsSquare != true)
            {
                throw new System.ArgumentException("Matrix A is not a square matrix.");
            }


            if (matrixA.RowCount != matrixB.Count)
            {
                throw new System.ArgumentException("Matrix dimensions are not valid.");
            }
        }
        private static void CheckDimensions(ComplexMatrix matrixA, ComplexVector matrixB)
        {
            if (matrixA.IsSquare != true)
            {
                throw new System.ArgumentException("Matrix A is not a square matrix.");
            }


            if (matrixA.RowCount != matrixB.Count)
            {
                throw new System.ArgumentException("Matrix dimensions are not valid.");
            }
        }

        private static void CheckDimensions(BaseMatrix matrixA, Vector vectorB)
        {
            if (!matrixA.IsSquare)
            {
                throw new System.ArgumentException("Matrix A is not a square matrix.");
            }

            if (vectorB.Type != VectorType.Column)
            {
                throw new System.ArgumentException("The vector should be a column vector.");
            }

            if (matrixA.RowCount != vectorB.Count)
            {
                throw new System.ArgumentException("Matrix dimensions are not valid.");
            }
        }


        #endregion



    }
}
