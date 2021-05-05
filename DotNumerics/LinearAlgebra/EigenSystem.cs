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
using DotNumerics.LinearAlgebra.CSLapack;
using DotNumerics.LinearAlgebra.CSEispack;
using DotNumerics.FortranLibrary;
using System.Linq;

namespace DotNumerics.LinearAlgebra
{
    /// <summary>
    /// Computes the eigenvalues and the eigenvectors of a square matrix.
    /// </summary>
    public sealed class EigenSystem
    {

        #region Fields

        DGEEV _dgeev;
        DSYEV _dsyev;
        DSBEV _dsbev;
        CG _cg;

        #endregion


        /// <summary>
        /// Initializes a new instance of the EigenSystem class.
        /// </summary>
        public EigenSystem()
        {

        }

        #region General matrix


        /// <summary>
        /// Computes the eigenvalues for an N-by-N real nonsymmetric matrix A.
        /// Modified By MRB, return real part of eigen values
        /// </summary>
        /// <param name="A">N-by-N real nonsymmetric matrix A.</param>
        /// <returns>The eigenvalues.</returns>
        public Matrix GetEigenvalues_Real(Matrix A)
        {

            if (this._dgeev == null) this._dgeev = new DGEEV();

            this.CheckDimensions(A);


            Matrix ACopy = A.Clone();
            double[] ACopyData = ACopy.Data;
            Matrix RealEVectors = new Matrix(1, 1);
            double[] EigenVectsData = RealEVectors.Data;
            Matrix EigenVals = new Matrix(A.RowCount, 1);

            double[] REigVal = new double[A.RowCount];
            double[] IEigVal = new double[A.RowCount];

            //double[] EigenValsData = EigenVals.Data;
            int Info = 0;

            double[] VL = new double[A.RowCount];

            double[] Work = new double[1];
            int LWork = -1;

            //Calculamos LWORK 
            _dgeev.Run("N", "N", A.RowCount, ref ACopyData, 0, ACopy.RowCount, ref REigVal, 0, ref IEigVal, 0, ref VL, 0, 1, ref EigenVectsData, 0, A.RowCount, ref Work, 0, LWork, ref Info);

            LWork = Convert.ToInt32(Work[0]);
            if (LWork > 0)
            {
                Work = new double[LWork];
                _dgeev.Run("N", "N", A.RowCount, ref ACopyData, 0, ACopy.RowCount, ref REigVal, 0, ref IEigVal, 0, ref VL, 0, 1, ref EigenVectsData, 0, A.RowCount, ref Work, 0, LWork, ref Info);
            }
            else
            {

                //Error
            }


            #region Error
            //= 0:  successful exit
            //.LT. 0:  if INFO = -i, the i-th argument had an illegal value.
            //.GT. 0:  if INFO = i, the QR algorithm failed to compute all the
            // eigenvalues, and no eigenvectors have been computed;
            // elements i+1:N of WR and WI contain eigenvalues which
            // have converged.

            if (Info < 0)
            {
                string infoSTg = Math.Abs(Info).ToString();
                throw new ArgumentException("the " + infoSTg + " -th argument had an illegal value");
            }
            else if (Info > 0)
            {
                string infoSTg = Math.Abs(Info).ToString();
                throw new Exception("The QR algorithm failed to compute all the eigenvalues.");
            }

            #endregion


            for (int i = 0; i < EigenVals.RowCount; i++)
            {
                EigenVals[i, 0] = REigVal[i];
            }


            return EigenVals;

        }

        /// <summary>
        /// Computes the eigenvalues for an N-by-N real nonsymmetric matrix A.
        /// </summary>
        /// <param name="A">N-by-N real nonsymmetric matrix A.</param>
        /// <returns>The eigenvalues.</returns>
        public ComplexVector GetEigenvalues(Matrix A)
        {

            if (this._dgeev == null) this._dgeev = new DGEEV();

            this.CheckDimensions(A);


            Matrix ACopy = A.Clone();
            double[] ACopyData = ACopy.Data;
            Matrix RealEVectors = new Matrix(1, 1);
            double[] EigenVectsData = RealEVectors.Data;
            ComplexVector EigenVals = new ComplexVector(A.RowCount);

            double[] REigVal = new double[A.RowCount];
            double[] IEigVal = new double[A.RowCount];

            //double[] EigenValsData = EigenVals.Data;
            int Info = 0;

            double[] VL = new double[A.RowCount];

            double[] Work = new double[1];
            int LWork = -1;

            //Calculamos LWORK 
            _dgeev.Run("N", "N", A.RowCount, ref ACopyData, 0, ACopy.RowCount, ref REigVal, 0, ref IEigVal, 0, ref VL, 0, 1, ref EigenVectsData, 0, A.RowCount, ref Work, 0, LWork, ref Info);

            LWork = Convert.ToInt32(Work[0]);
            if (LWork > 0)
            {
                Work = new double[LWork];
                _dgeev.Run("N", "N", A.RowCount, ref ACopyData, 0, ACopy.RowCount, ref REigVal, 0, ref IEigVal, 0, ref VL, 0, 1, ref EigenVectsData, 0, A.RowCount, ref Work, 0, LWork, ref Info);
            }
            else
            {

                //Error
            }


            #region Error
            //= 0:  successful exit
            //.LT. 0:  if INFO = -i, the i-th argument had an illegal value.
            //.GT. 0:  if INFO = i, the QR algorithm failed to compute all the
            // eigenvalues, and no eigenvectors have been computed;
            // elements i+1:N of WR and WI contain eigenvalues which
            // have converged.

            if (Info < 0)
            {
                string infoSTg = Math.Abs(Info).ToString();
                throw new ArgumentException("the " + infoSTg + " -th argument had an illegal value");
            }
            else if (Info > 0)
            {
                string infoSTg = Math.Abs(Info).ToString();
                throw new Exception("The QR algorithm failed to compute all the eigenvalues.");
            }

            #endregion


            for (int i = 0; i < EigenVals.Count; i++)
            {
                EigenVals[i] = new Complex(REigVal[i], IEigVal[i]);
            }


            return EigenVals;

        }

        /// <summary>
        /// Computes for an N-by-N real nonsymmetric matrix A, the
        /// eigenvalues and eigenvectors.
        /// Modified By MRB, return real part of eigen values and eigen vectors
        /// </summary>
        /// <param name="A">N-by-N real nonsymmetric matrix A.</param>
        /// <param name="EigenVectors">The eigenvectors.</param>
        /// <returns>The eigenvalues as a column (Nx1)</returns>
        public Matrix GetEigenvalues_Real(Matrix A, out Matrix EigenVectors)
        {

            if (this._dgeev == null) this._dgeev = new DGEEV();

            this.CheckDimensions(A);


            Matrix ACopy = A.Clone();
            double[] ACopyData = ACopy.Data;
            EigenVectors = new Matrix(A.RowCount, A.ColumnCount);
            Matrix RealEVectors = new Matrix(A.RowCount, A.ColumnCount);
            double[] EigenVectsData = RealEVectors.Data;
            Matrix EigenVals = new Matrix(A.RowCount, 1);

            double[] REigVal = new double[A.RowCount];
            double[] IEigVal = new double[A.RowCount];

            //double[] EigenValsData = EigenVals.Data;
            int Info = 0;

            double[] VL = new double[A.RowCount];

            double[] Work = new double[1];
            int LWork = -1;

            //Calculamos LWORK 
            _dgeev.Run("N", "V", A.RowCount, ref ACopyData, 0, ACopy.RowCount, ref REigVal, 0, ref IEigVal, 0, ref VL, 0, 1, ref EigenVectsData, 0, A.RowCount, ref Work, 0, LWork, ref Info);

            LWork = Convert.ToInt32(Work[0]);
            if (LWork > 0)
            {
                Work = new double[LWork];
                _dgeev.Run("N", "V", A.RowCount, ref ACopyData, 0, ACopy.RowCount, ref REigVal, 0, ref IEigVal, 0, ref VL, 0, 1, ref EigenVectsData, 0, A.RowCount, ref Work, 0, LWork, ref Info);
            }
            else
            {

                //Error
            }


            #region Error
            //= 0:  successful exit
            //.LT. 0:  if INFO = -i, the i-th argument had an illegal value.
            //.GT. 0:  if INFO = i, the QR algorithm failed to compute all the
            // eigenvalues, and no eigenvectors have been computed;
            // elements i+1:N of WR and WI contain eigenvalues which
            // have converged.

            if (Info < 0)
            {
                string infoSTg = Math.Abs(Info).ToString();
                throw new ArgumentException("the " + infoSTg + " -th argument had an illegal value");
            }
            else if (Info > 0)
            {
                string infoSTg = Math.Abs(Info).ToString();
                throw new Exception("The QR algorithm failed to compute all the eigenvalues.");
            }

            #endregion


            for (int i = 0; i < EigenVals.RowCount; i++)
            {
                EigenVals[i, 0] = REigVal[i];
                if (IEigVal[i] >= 0)
                    for (int j = 0; j < EigenVectors.RowCount; j++)
                    {
                        EigenVectors[j, i] = RealEVectors[j, i];
                    }
                else
                    for (int j = 0; j < EigenVectors.RowCount; j++)
                    {
                        EigenVectors[j, i] = RealEVectors[j, i - 1];
                    }
            }

            return EigenVals;
        }

        /// <summary>
        /// Computes for an N-by-N real nonsymmetric matrix A, the
        /// eigenvalues and eigenvectors.
        /// </summary>
        /// <param name="A">N-by-N real nonsymmetric matrix A.</param>
        /// <param name="EigenVectors">The eigenvectors.</param>
        /// <returns>The eigenvalues.</returns>
        public ComplexMatrix GetEigenvalues(Matrix A, out ComplexMatrix EigenVectors)
        {

            if (this._dgeev == null) this._dgeev = new DGEEV();

            this.CheckDimensions(A);


            Matrix ACopy = A.Clone();
            double[] ACopyData = ACopy.Data;
            EigenVectors = new ComplexMatrix(A.RowCount, A.ColumnCount);
            Matrix RealEVectors = new Matrix(A.RowCount, A.ColumnCount);
            double[] EigenVectsData = RealEVectors.Data;
            ComplexMatrix EigenVals = new ComplexMatrix(A.RowCount, 1);

            double[] REigVal = new double[A.RowCount];
            double[] IEigVal = new double[A.RowCount];

            //double[] EigenValsData = EigenVals.Data;
            int Info = 0;

            double[] VL = new double[A.RowCount];

            double[] Work = new double[1];
            int LWork = -1;

            //Calculamos LWORK 
            _dgeev.Run("N", "V", A.RowCount, ref ACopyData, 0, ACopy.RowCount, ref REigVal, 0, ref IEigVal, 0, ref VL, 0, 1, ref EigenVectsData, 0, A.RowCount, ref Work, 0, LWork, ref Info);

            LWork = Convert.ToInt32(Work[0]);
            if (LWork > 0)
            {
                Work = new double[LWork];
                _dgeev.Run("N", "V", A.RowCount, ref ACopyData, 0, ACopy.RowCount, ref REigVal, 0, ref IEigVal, 0, ref VL, 0, 1, ref EigenVectsData, 0, A.RowCount, ref Work, 0, LWork, ref Info);
            }
            else
            {

                //Error
            }


            #region Error
            //= 0:  successful exit
            //.LT. 0:  if INFO = -i, the i-th argument had an illegal value.
            //.GT. 0:  if INFO = i, the QR algorithm failed to compute all the
            // eigenvalues, and no eigenvectors have been computed;
            // elements i+1:N of WR and WI contain eigenvalues which
            // have converged.

            if (Info < 0)
            {
                string infoSTg = Math.Abs(Info).ToString();
                throw new ArgumentException("the " + infoSTg + " -th argument had an illegal value");
            }
            else if (Info > 0)
            {
                string infoSTg = Math.Abs(Info).ToString();
                throw new Exception("The QR algorithm failed to compute all the eigenvalues.");
            }

            #endregion


            for (int i = 0; i < EigenVals.RowCount; i++)
            {
                EigenVals[i, 0] = new Complex(REigVal[i], IEigVal[i]);
            }

            for (int i = 0; i < EigenVals.RowCount; i++)
            {
                if (EigenVals[i, 0].Imaginary == 0.0)
                {
                    for (int j = 0; j < EigenVectors.RowCount; j++)
                    {
                        EigenVectors[j, i] = new Complex(RealEVectors[j, i], 0.0);
                    }
                }
                else
                {
                    if (EigenVals[i, 0].Imaginary > 0.0)
                    {
                        for (int j = 0; j < EigenVectors.RowCount; j++)
                        {
                            EigenVectors[j, i] = new Complex(RealEVectors[j, i], RealEVectors[j, i + 1]);
                        }
                    }
                    else
                    {
                        for (int j = 0; j < EigenVectors.RowCount; j++)
                        {
                            EigenVectors[j, i] = new Complex(RealEVectors[j, i - 1], -RealEVectors[j, i]);

                        }
                    }
                }
            }

            return EigenVals;
        }

        #endregion


        #region SymmetricMatrix

        /// <summary>
        /// Computes all eigenvalues of a real symmetric matrix A.
        /// </summary>
        /// <param name="A">The real symmetric matrix A.</param>
        /// <returns>The eigenvalues.</returns>
        public Matrix GetEigenvalues(SymmetricMatrix A)
        {
            if (this._dsyev == null) this._dsyev = new DSYEV();


            this.CheckDimensions(A);

            Matrix EigenVects = new Matrix(A.RowCount, A.ColumnCount, A.Data);
            double[] EigenVectsData = EigenVects.Data;
            Matrix EigenVals = new Matrix(A.RowCount, 1);
            double[] EigenValsData = EigenVals.Data;
            int Info = 0;

            double[] Work = new double[1];

            int LWork = -1;
            //Calculamos LWORK ideal

            _dsyev.Run("N", "U", A.RowCount, ref EigenVectsData, 0, A.RowCount, ref EigenValsData, 0, ref Work, 0, LWork, ref Info);

            LWork = Convert.ToInt32(Work[0]);
            if (LWork > 0)
            {
                Work = new double[LWork];
                _dsyev.Run("N", "U", A.RowCount, ref EigenVectsData, 0, A.RowCount, ref EigenValsData, 0, ref Work, 0, LWork, ref Info);
            }
            else
            {

                //Error
            }

            #region Error
            /// = 0:  successful exit
            /// .LT. 0:  if INFO = -i, the i-th argument had an illegal value
            /// .GT. 0:  if INFO = i, the algorithm failed to converge; i
            /// off-diagonal elements of an intermediate tridiagonal
            /// form did not converge to zero.

            if (Info < 0)
            {
                string infoSTg = Math.Abs(Info).ToString();
                throw new ArgumentException("the " + infoSTg + " -th argument had an illegal value");
            }
            else if (Info > 0)
            {
                string infoSTg = Math.Abs(Info).ToString();
                throw new Exception("The algorithm failed to converge.");
            }

            #endregion


            return EigenVals;
        }

        /// <summary>
        /// Computes all eigenvalues and eigenvectors of a of a real symmetric matrix A.
        /// </summary>
        /// <param name="A">The real symmetric matrix A.</param>
        /// <param name="EigenVects">The eigenvectors.</param>
        /// <returns>The eigenvalues.</returns>
        public Matrix GetEigenvalues(SymmetricMatrix A, out Matrix EigenVects)
        {

            if (this._dsyev == null) this._dsyev = new DSYEV();
            this.CheckDimensions(A);


            EigenVects = new Matrix(A.RowCount, A.ColumnCount, A.Data);
            double[] EigenVectsData = EigenVects.Data;
            Matrix EigenVals = new Matrix(A.RowCount, 1);
            double[] EigenValsData = EigenVals.Data;
            int Info = 0;

            double[] Work = new double[1];

            int LWork = -1;
            //Calculamos LWORK ideal

            _dsyev.Run("V", "U", A.RowCount, ref EigenVectsData, 0, A.RowCount, ref EigenValsData, 0, ref Work, 0, LWork, ref Info);

            LWork = Convert.ToInt32(Work[0]);
            if (LWork > 0)
            {
                Work = new double[LWork];
                _dsyev.Run("V", "U", A.RowCount, ref EigenVectsData, 0, A.RowCount, ref EigenValsData, 0, ref Work, 0, LWork, ref Info);
            }
            else
            {

                //Error
            }

            #region Error
            /// = 0:  successful exit
            /// .LT. 0:  if INFO = -i, the i-th argument had an illegal value
            /// .GT. 0:  if INFO = i, the algorithm failed to converge; i
            /// off-diagonal elements of an intermediate tridiagonal
            /// form did not converge to zero.

            if (Info < 0)
            {
                string infoSTg = Math.Abs(Info).ToString();
                throw new ArgumentException("the " + infoSTg + " -th argument had an illegal value");
            }
            else if (Info > 0)
            {
                string infoSTg = Math.Abs(Info).ToString();
                throw new Exception("The algorithm failed to converge.");
            }

            #endregion


            return EigenVals;
        }

        #endregion


        #region SymmetricBandMatrix

        /// <summary>
        ///Computes all the eigenvalues of
        /// a real symmetric band matrix A.
        /// </summary>
        /// <param name="A">The real symmetric band matrix A.</param>
        /// <returns>The eigenvalues.</returns>
        public Matrix GetEigenvalues(SymmetricBandMatrix A)
        {
            if (this._dsbev == null) this._dsbev = new DSBEV();

            this.CheckDimensions(A);

            Matrix SymmetricBand = A.GetSymmetricBandPackedMatrix();
            double[] SymmetricBandData = SymmetricBand.Data;
            Matrix EigenVects = new Matrix(1, 1);  //Se pone (1,1) pues nos e usaran
            double[] EigenVectsData = EigenVects.Data;
            Matrix EigenVals = new Matrix(A.RowCount, 1);
            double[] EigenValsData = EigenVals.Data;
            int Info = 0;

            double[] Work = new double[3 * A.RowCount - 2];


            _dsbev.Run("N", "U", A.RowCount, A.UpperBandWidth, ref SymmetricBandData, 0, SymmetricBand.RowCount, ref EigenValsData, 0, ref EigenVectsData, 0, A.RowCount, ref Work, 0, ref Info);


            #region Error
            /// = 0:  successful exit
            /// .LT. 0:  if INFO = -i, the i-th argument had an illegal value
            /// .GT. 0:  if INFO = i, the algorithm failed to converge; i
            /// off-diagonal elements of an intermediate tridiagonal
            /// form did not converge to zero.

            if (Info < 0)
            {
                string infoSTg = Math.Abs(Info).ToString();
                throw new ArgumentException("the " + infoSTg + " -th argument had an illegal value");
            }
            else if (Info > 0)
            {
                string infoSTg = Math.Abs(Info).ToString();
                throw new Exception("The algorithm failed to converge.");
            }

            #endregion


            return EigenVals;
        }

        /// <summary>
        ///Computes all the eigenvalues and eigenvectors of
        /// a real symmetric band matrix A.
        /// </summary>
        /// <param name="A">The real symmetric band matrix A.</param>
        /// <param name="EigenVects">The eigenvectors.</param>
        /// <returns>The eigenvalues.</returns>
        public Matrix GetEigenvalues(SymmetricBandMatrix A, out Matrix EigenVects)
        {
            if (this._dsbev == null) this._dsbev = new DSBEV();
            this.CheckDimensions(A);

            Matrix SymmetricBand = A.GetSymmetricBandPackedMatrix();
            double[] SymmetricBandData = SymmetricBand.Data;
            EigenVects = new Matrix(A.RowCount, A.ColumnCount);
            double[] EigenVectsData = EigenVects.Data;
            Matrix EigenVals = new Matrix(A.RowCount, 1);
            double[] EigenValsData = EigenVals.Data;
            int Info = 0;

            double[] Work = new double[3 * A.RowCount - 2];


            _dsbev.Run("V", "U", A.RowCount, A.UpperBandWidth, ref SymmetricBandData, 0, SymmetricBand.RowCount, ref EigenValsData, 0, ref EigenVectsData, 0, A.RowCount, ref Work, 0, ref Info);


            #region Error
            /// = 0:  successful exit
            /// .LT. 0:  if INFO = -i, the i-th argument had an illegal value
            /// .GT. 0:  if INFO = i, the algorithm failed to converge; i
            /// off-diagonal elements of an intermediate tridiagonal
            /// form did not converge to zero.

            if (Info < 0)
            {
                string infoSTg = Math.Abs(Info).ToString();
                throw new ArgumentException("the " + infoSTg + " -th argument had an illegal value");
            }
            else if (Info > 0)
            {
                string infoSTg = Math.Abs(Info).ToString();
                throw new Exception("The algorithm failed to converge.");
            }

            #endregion


            return EigenVals;
        }

        #endregion


        #region Complex General matrix

        /// <summary>
        /// Computes the eigenvalues for an complex general matrix A.
        /// </summary>
        /// <param name="A">The complex general matrix A.</param>
        /// <returns>The eigenvalues.</returns>
        public ComplexVector GetEigenvalues(ComplexMatrix A)
        {
            //Fortran Ejemplo
            //CG(NM,N,AR,AI,WR,WI,1,ZR,ZI,SCALE,ORTR,ORTI,ERROR)
            //C
            //C     THIS DRIVER TESTS  EISPACK  FOR THE CLASS OF COMPLEX GENERAL
            //C     MATRICES SUMMARIZING THE FIGURES OF MERIT FOR ALL PATHS.
            //C
            //C     THIS DRIVER IS CATALOGUED AS  EISPDRV4(CGSUMARY).
            //C
            //C     THE DIMENSION OF  AR,AI,ZR,ZI,ASAVER,ASAVEI,RM1,  AND  RM2 SHOULD
            //C     BE  NM  BY  NM.
            //C     THE DIMENSION OF  WR,WI,WR1,WI1,SELECT,SLHOLD,INT,SCALE,ORTR,ORTI,
            //C     RV1  AND  RV2  SHOULD BE  NM.
            //C     THE DIMENSION OF  ARHOLD  AND  AIHOLD  SHOULD BE  NM  BY  NM.
            //C     HERE NM = 20.


            if (this._cg == null) this._cg = new CG();

            this.CheckDimensions(A);

            Matrix AReal = A.GetReal();
            double[] ARealData = AReal.Data;
            Matrix AImag = A.GetImag();
            double[] AImagData = AImag.Data;

            ComplexMatrix EigenVals = new ComplexMatrix(A.RowCount, 1);
            Matrix RealEigenVals = new Matrix(A.RowCount, 1);
            double[] RealEigenValsData = RealEigenVals.Data;
            Matrix ImagEigenVals = new Matrix(A.RowCount, 1);
            double[] ImagEigenValsData = ImagEigenVals.Data;

            ComplexMatrix EigenVectors = new ComplexMatrix(1, 1);
            Matrix RealEigVect = new Matrix(A.RowCount);
            double[] RealEigVectData = RealEigVect.Data;
            Matrix ImagEigVect = new Matrix(A.RowCount);
            double[] ImagEigVectData = ImagEigVect.Data;

            double[] SCALE = new double[A.RowCount];
            double[] ORTR = new double[A.RowCount];
            double[] ORTI = new double[A.RowCount];

            int Info = 0;

            int matz = 0; //No se calculan los eigenvectores

            _cg.Run(A.RowCount, A.RowCount, ref ARealData, 0, ref AImagData, 0, ref RealEigenValsData, 0, ref ImagEigenValsData, 0,
                matz, ref RealEigVectData, 0, ref ImagEigVectData, 0, ref SCALE, 0, ref ORTR, 0, ref ORTI, 0, ref Info);


            #region Error
            /// is set to
            /// zero       for normal return,
            /// j          if the limit of 30*n iterations is exhausted
            /// while the j-th eigenvalue is being sought.

            if (Info != 0)
            {
                throw new ArgumentException("The limit of 30*n iterations is exhausted");
            }

            #endregion



            EigenVals.SetReal(RealEigenVals);
            EigenVals.SetImag(ImagEigenVals);

            //EigenVectors.SetReal(RealEigVect);
            //EigenVectors.SetImag(ImagEigVect);

            return EigenVals.GetColumnVectors()[0];
        }

        /// <summary>
        /// Computes the eigenvalues and eigenvectors for an complex general matrix A.
        /// </summary>
        /// <param name="A">The complex general matrix A.</param>
        /// <param name="EigenVectors">The eigenvectors.</param>
        /// <returns>The eigenvalues.</returns>
        public ComplexVector GetEigenvalues(ComplexMatrix A, out ComplexMatrix EigenVectors)
        {
            //Fortran Ejemplo
            //CG(NM,N,AR,AI,WR,WI,1,ZR,ZI,SCALE,ORTR,ORTI,ERROR)
            //C
            //C     THIS DRIVER TESTS  EISPACK  FOR THE CLASS OF COMPLEX GENERAL
            //C     MATRICES SUMMARIZING THE FIGURES OF MERIT FOR ALL PATHS.
            //C
            //C     THIS DRIVER IS CATALOGUED AS  EISPDRV4(CGSUMARY).
            //C
            //C     THE DIMENSION OF  AR,AI,ZR,ZI,ASAVER,ASAVEI,RM1,  AND  RM2 SHOULD
            //C     BE  NM  BY  NM.
            //C     THE DIMENSION OF  WR,WI,WR1,WI1,SELECT,SLHOLD,INT,SCALE,ORTR,ORTI,
            //C     RV1  AND  RV2  SHOULD BE  NM.
            //C     THE DIMENSION OF  ARHOLD  AND  AIHOLD  SHOULD BE  NM  BY  NM.
            //C     HERE NM = 20.

            if (this._cg == null) this._cg = new CG();

            this.CheckDimensions(A);

            Matrix AReal = A.GetReal();
            double[] ARealData = AReal.Data;
            Matrix AImag = A.GetImag();
            double[] AImagData = AImag.Data;

            ComplexMatrix EigenVals = new ComplexMatrix(A.RowCount, 1);
            Matrix RealEigenVals = new Matrix(A.RowCount, 1);
            double[] RealEigenValsData = RealEigenVals.Data;
            Matrix ImagEigenVals = new Matrix(A.RowCount, 1);
            double[] ImagEigenValsData = ImagEigenVals.Data;

            EigenVectors = new ComplexMatrix(A.RowCount, A.ColumnCount);
            Matrix RealEigVect = new Matrix(A.RowCount);
            double[] RealEigVectData = RealEigVect.Data;
            Matrix ImagEigVect = new Matrix(A.RowCount);
            double[] ImagEigVectData = ImagEigVect.Data;

            double[] SCALE = new double[A.RowCount];
            double[] ORTR = new double[A.RowCount];
            double[] ORTI = new double[A.RowCount];

            int Info = 0;
            int matz = 1; //Se calculan los eigenvalores y los eigenvectores
            _cg.Run(A.RowCount, A.RowCount, ref ARealData, 0, ref AImagData, 0, ref RealEigenValsData, 0, ref ImagEigenValsData, 0,
                matz, ref RealEigVectData, 0, ref ImagEigVectData, 0, ref SCALE, 0, ref ORTR, 0, ref ORTI, 0, ref Info);


            #region Error
            /// is set to
            /// zero       for normal return,
            /// j          if the limit of 30*n iterations is exhausted
            /// while the j-th eigenvalue is being sought.

            if (Info != 0)
            {
                throw new ArgumentException("The limit of 30*n iterations is exhausted");
            }

            #endregion



            EigenVals.SetReal(RealEigenVals);
            EigenVals.SetImag(ImagEigenVals);

            EigenVectors.SetReal(RealEigVect);
            EigenVectors.SetImag(ImagEigVect);

            return EigenVals.GetColumnVectors()[0];
        }



        #endregion


        #region Private methods

        private void CheckDimensions(BaseMatrix matrixA)
        {
            if (matrixA.IsSquare != true)
            {
                throw new System.ArgumentException("Matrix A is not a square matrix.");
            }

        }

        private void CheckDimensions(ComplexMatrix matrixA)
        {
            if (matrixA.IsSquare != true)
            {
                throw new System.ArgumentException("Matrix A is not a square matrix.");
            }

        }

        #endregion

        #region Ease of access, By MRB

        #region tools
        public static void Sort(Vector Val, Matrix Vec, int num, bool ascending)
        {
            for (int i = 0; i < num; i++)
            {
                for (int j = i + 1; j < Val.Count; j++)
                    if (ascending)
                    {
                        if (Val[i] > Val[j])
                        {
                            Val.Swap(i, j);
                            Vec.SwapColumns(i, j);
                        }
                    }
                    else
                    {
                        if (Val[i] < Val[j])
                        {
                            Val.Swap(i, j);
                            Vec.SwapColumns(i, j);
                        }
                    }
            }
        }
        public enum SortMode
        {
            Modulus, Phase, Real, Imag, Modulus_then_neg_imag, FrequencyInStateSpace
        }
        public static void Sort(ComplexVector Val, ComplexMatrix Vec, int num, bool ascending, SortMode sortMode)
        {
            for (int i = 0; i < num; i++)
            {
                for (int j = i + 1; j < Val.Count; j++)
                    if (MustSwap(Val[i], Val[j], ascending, sortMode))
                    {
                        Val.Swap(i, j);
                        Vec.SwapColumns(i, j);
                    }
            }
        }
        static bool MustSwap(Complex a, Complex b, bool ascending, SortMode sortMode)
        {
            if (sortMode == SortMode.Modulus_then_neg_imag)
            {
                if (ascending)
                {
                    if (a.Modulus > b.Modulus)
                        if (b.Imaginary < 0)
                            return true;
                }
                else
                {
                    if (a.Modulus < b.Modulus)
                        if (b.Imaginary < 0)
                            return true;
                }
            }
            else if (sortMode == SortMode.Modulus)
            {
                var a_ = a.Modulus;
                var b_ = b.Modulus;
                if (ascending)
                {
                    if (a_ > b_)
                        return true;
                    if (Math.Abs(a_ - b_) / (a_ + 1e-100) < 1e-8)
                        return a.Imaginary + a.Real > b.Imaginary + b.Real;

                }
                else
                {
                    if (a.Modulus < b.Modulus)
                        return true;
                    if (Math.Abs(a_ - b_) / (a_ + 1e-100) < 1e-8)
                        return a.Imaginary + a.Real < b.Imaginary + b.Real;
                }
            }
            else if (sortMode == SortMode.Imag)
            {
                if (ascending)
                {
                    if (Math.Abs(a.Imaginary - b.Imaginary) / (Math.Abs(a.Imaginary) + 1e-100) < 1e-8)
                        return a.Real > b.Real;
                    else
                        return a.Imaginary > b.Imaginary;
                }
                else
                {
                    if (a.Imaginary < b.Imaginary)
                        return true;
                }
            }
            else if (sortMode == SortMode.Real)
            {
                if (ascending)
                {
                    if (a.Real > b.Real)
                        return true;
                }
                else
                {
                    if (a.Real < b.Real)
                        return true;
                }
            }
            else if (sortMode == SortMode.Phase)
            {
                if (ascending)
                {
                    if (a.Angle > a.Angle)
                        return true;
                }
                else
                {
                    if (a.Angle < b.Angle)
                        return true;
                }
            }
            return false;
        }
        public static void Sort(Vector Val, int num, bool ascending)
        {
            for (int i = 0; i < num; i++)
            {
                for (int j = i + 1; j < Val.Count; j++)
                    if (ascending)
                    {
                        if (Val[i] > Val[j])
                            Val.Swap(i, j);
                    }
                    else
                        if (Val[i] < Val[j])
                        Val.Swap(i, j);
            }
        }
        public static void Sort(ComplexVector Val, int num, bool ascending, SortMode sortMode)
        {
            for (int i = 0; i < num; i++)
                for (int j = i + 1; j < Val.Count; j++)
                    if (MustSwap(Val[i], Val[j], ascending, sortMode))
                        Val.Swap(i, j);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="lambda">sigma +- i omega</param>
        /// <param name="ModeShapes">هر ستون یک شکل مود</param>
        public static void SortFrequencyStateSpace(ref ComplexVector lambda, ref ComplexMatrix ModeShapes, int num)
        {
            var real_count = 0;
            for (int i = 0; i < lambda.Count; i++)
                if (lambda[i].Imaginary == 0 && lambda[i].Real > 1e-9)
                    real_count++;
            var n_ = lambda.Count / 2 - real_count / 2;
            if (num < 0) num = n_;
            Sort(lambda, ModeShapes, lambda.Count, true, SortMode.Imag);
            /*var N = 0;
            for (; N < num / 2; N++)
                if (Math.Abs(lambda[n_ + N].Imaginary) + Math.Abs(lambda[n_ - N - 1].Imaginary) > 1e-10)
                    break;
            var j = n_ - N;
            for (int i = n_; i < n_ + N; i++)
            {
                lambda.Swap(i, j);
                ModeShapes.SwapColumns(i, j);
                j += 2;
            }*/
            lambda = lambda[n_, num];
            ModeShapes = ModeShapes[0, n_, -1, num];
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="lambda">sigma +- i omega</param> 
        public static void SortFrequencyStateSpace(ref ComplexVector lambda, int num)
        {
            var n_ = lambda.Count / 2;
            if (num < 0) num = n_;
            Sort(lambda, lambda.Count, true, SortMode.Imag);
            /*var N = 0;
            for (; N < n_; N++)
                if (Math.Abs(lambda[n_ + N].Imaginary) + Math.Abs(lambda[n_ - N - 1].Imaginary) > 1e-10)
                    break;
            var j = n_ - N;
            for (int i = n_; i < n_ + N; i++)
            {
                lambda.Swap(i, j);
                j += 2;
            }*/
            lambda = lambda[n_, num];
        }
        public static void Normalize(Matrix Vec, Matrix M, bool massNormalize, bool normalize)
        {
            var num = Vec.ColumnCount;
            if (normalize)
                massNormalize = false;
            if (massNormalize)
                for (int i = 0; i < num; i++)
                {
                    var V = Vec.GetColumnVector(i);
                    var r = V.DotProduct(M * V);
                    if (r == 0) r = 1;
                    V /= Math.Sqrt(r);
                    Vec.InsertAt(0, i, V, -1, 1);
                }
            if (normalize)
                for (int i = 0; i < num; i++)
                {
                    var V = Vec.GetColumnVector(i);
                    V.NormalizeInPlace();
                    Vec.InsertAt(0, i, V, -1, 1);
                }
        }
        public static void Normalize(ComplexMatrix Vec, Matrix M, bool massNormalize, bool normalize)
        {
            var num = Vec.ColumnCount;
            if (normalize)
                massNormalize = false;
            if (massNormalize)
                for (int i = 0; i < num; i++)
                {
                    var V = Vec.SubMatrix(startCol: i, colCount: 1);
                    var r = (V.Transpose() * (M * V))[0, 0];
                    if (r == 0) r = 1;
                    V /= Math.Sqrt(r.Modulus);
                    Vec.InsertAt(0, i, V, -1, 1);
                }
            if (normalize)
                for (int i = 0; i < num; i++)
                {
                    var V = Vec.GetColumnVector(i);
                    V /= V.Norm();
                    Vec.InsertAt(0, i, V, -1, 1);
                }
        }
        public static void Normalize(ComplexMatrix Vec, ComplexMatrix M, bool massNormalize, bool normalize)
        {
            var num = Vec.ColumnCount;
            if (normalize)
                massNormalize = false;
            if (massNormalize)
                for (int i = 0; i < num; i++)
                {
                    var V = Vec.SubMatrix(startCol: i, colCount: 1);
                    var r = (V.Transpose() * (M * V))[0, 0];
                    if (r == 0) r = 1;
                    V /= Math.Sqrt(r.Modulus);
                    Vec.InsertAt(0, i, V, -1, 1);
                }
            if (normalize)
                for (int i = 0; i < num; i++)
                {
                    var V = Vec.GetColumnVector(i);
                    V /= V.Norm();
                    Vec.InsertAt(0, i, V, -1, 1);
                }
        }
        #endregion
        /// <summary>
        /// Added by MRB
        /// K X = λ M X
        /// </summary>
        /// <param name="K"></param>
        /// <param name="M"></param>
        /// <param name="Vec"></param>
        /// <param name="sort"></param>
        /// <param name="num"></param>
        /// <returns></returns>
        public static Vector Eig(Matrix K, Matrix M, out Matrix Vec, bool sort = true, int num = -1, bool massNormalize = true, bool normalize = false, bool ascending = true)
        {
            if (num <= 0)
                num = M.ColumnCount;
            Vector Val = null;
            if (ARMADILO.Enabled)
            {
                ComplexMatrix Vec_;
                var Val_ = ARMADILO.eig_pair(K, M, out Vec_);
                Val = Val_.Real();
                Vec = Vec_.Real();
            }
            else
                Val = (new EigenSystem()).GetEigenvalues_Real(M.Inverse() * K, out Vec).GetColumnVector(0);
            if (sort)
                Sort(Val, Vec, num, ascending);
            if (ExportDebugData)
            {
                Vec.Export4Matlab("vec.m");
                Val.Export4Matlab("val.m");
            }
            Vec = Vec.SubMatrix(0, 0, colCount: num);
            Val.Resize(num);
            Normalize(Vec, M, massNormalize, normalize);
            return Val;
        }
        /// <summary>
        /// Added by MRB
        /// </summary>
        /// <param name="K"></param>
        /// <param name="M"></param>
        /// <param name="Vec"></param>
        /// <param name="sort"></param>
        /// <param name="num"></param>
        /// <returns></returns>
        public static ComplexVector Eig_(Matrix K, Matrix M, out ComplexMatrix Vec, bool sort = true, int num = -1, bool massNormalize = true, bool normalize = false, bool ascending = true, SortMode sortMode = SortMode.Modulus)
        {
            if (num <= 0)
                num = M.ColumnCount;
            ComplexVector Val;
            if (ARMADILO.Enabled)
                Val = ARMADILO.eig_pair(K, M, out Vec);
            else
                Val = (new EigenSystem()).GetEigenvalues(M.Inverse() * K, out Vec).GetColumnVector(0);
            if (sort)
                Sort(Val, Vec, num, ascending, sortMode);
            if (ExportDebugData)
            {
                Vec.Export4Matlab("vec.m");
                Val.Export4Matlab("val.m");
            }
            Vec = Vec.SubMatrix(0, 0, colCount: num);
            Val.Resize(num);
            Normalize(Vec, M, massNormalize, normalize);
            return Val;
        }
        /// <summary>
        /// Added by MRB
        /// </summary>
        /// <param name="K"></param>
        /// <param name="M"></param> 
        /// <param name="sort"></param>
        /// <param name="num"></param>
        /// <returns></returns>
        public static Vector Eig(Matrix K, Matrix M, bool sort = true, int num = -1, bool ascending = true)
        {
            if (num <= 0)
                num = M.ColumnCount;
            Matrix Vec0;
            Vector Val;
            if (ARMADILO.Enabled)
            {
                var Val_ = ARMADILO.eig_pair(K, M, out _);
                Val = Val_.Real();
            }
            else
                Val = (new EigenSystem()).GetEigenvalues_Real(M.Inverse() * K, out Vec0).GetColumnVector(0);
            if (sort)
                Sort(Val, num, ascending);
            if (ExportDebugData)
            {
                //Vec.Export4Matlab("vec.m");
                Val.Export4Matlab("val.m");
            }
            Val.Resize(num);
            return Val;
        }
        public static bool ExportDebugData = false;
        /// <summary>
        /// Added by MRB
        /// </summary>
        /// <param name="K"></param>
        /// <param name="M"></param> 
        /// <param name="sort"></param>
        /// <param name="num"></param>
        /// <returns></returns>
        public static ComplexVector Eig_(Matrix K, Matrix M, bool sort = true, int num = -1, bool ascending = true, SortMode sortMode = SortMode.Modulus)
        {
            if (num <= 0)
                num = M.ColumnCount;
            ComplexVector Val;
            if (ARMADILO.Enabled)
                Val = ARMADILO.eig_pair(K, M, out _);
            else
                Val = (new EigenSystem()).GetEigenvalues(M.Inverse() * K, out _).GetColumnVector(0);
            if (sort)
                Sort(Val, num, ascending, sortMode);
            if (ExportDebugData)
            {
                //Vec.Export4Matlab("vec.m");
                Val.Export4Matlab("val.m");
            }
            Val.Resize(num);
            return Val;
        }
        /// <summary>
        /// q = [ x_dot , x] ; 
        /// [M 0; 0 I] q_dot + [C K; -I 0] = 0*[f; 0] q
        /// </summary>
        /// <param name="K"></param>
        /// <param name="M"></param>
        /// <param name="C"></param>
        /// <param name="sort"></param>
        /// <param name="num"></param>
        /// <param name="ascending"></param>
        /// <param name="sortMode"></param>
        /// <returns></returns>
        public static ComplexVector Eig(Matrix K, Matrix M, Matrix C,
            bool sort = true, int num = -1, bool ascending = true,
            SortMode sortMode = SortMode.FrequencyInStateSpace)
        {
            var n = K.RowCount;
            var KK = new Matrix(2 * n, 2 * n);
            var MM = new Matrix(2 * n, 2 * n);
            var In = Matrix.I(n);
            KK.InsertAt(0, 0, C);
            KK.InsertAt(0, n, K);
            KK.InsertAt(n, 0, -In);
            MM.InsertAt(0, 0, M);
            MM.InsertAt(n, n, In);
            if (sort && sortMode == SortMode.FrequencyInStateSpace)
            {
                var lamb = Eig_(KK, MM, false, 2 * n);
                SortFrequencyStateSpace(ref lamb, num);
                return lamb;
            }
            return Eig_(KK, MM, sort, num, ascending, sortMode);
        }
        /// <summary>
        /// q = [ x_dot , x]  ;
        /// [M 0; 0 I] q_dot + [C K; -I 0] = 0*[f; 0] q
        /// </summary>
        /// <param name="K"></param>
        /// <param name="M"></param>
        /// <param name="C"></param>
        /// <param name="Vec"></param>
        /// <param name="sort"></param>
        /// <param name="num"></param>
        /// <param name="massNormalize"></param>
        /// <param name="normalize"></param>
        /// <param name="ascending"></param>
        /// <param name="sortMode"></param>
        /// <returns></returns>
        public static ComplexVector Eig(Matrix K, Matrix M, Matrix C,
            out ComplexMatrix Vec, bool sort = true,
            int num = -1, bool massNormalize = true,
            bool normalize = false, bool ascending = true,
            SortMode sortMode = SortMode.FrequencyInStateSpace)
        {
            var n = K.RowCount;
            var KK = new Matrix(2 * n, 2 * n);
            var MM = new Matrix(2 * n, 2 * n);
            var In = Matrix.I(n);
            KK.InsertAt(0, 0, C);
            KK.InsertAt(0, n, K);
            KK.InsertAt(n, 0, -In);
            MM.InsertAt(0, 0, M);
            MM.InsertAt(n, n, In);
            if (sort && sortMode == SortMode.FrequencyInStateSpace)
            {
                var lamb = Eig_(KK, MM, out Vec, false, 2 * n);
                SortFrequencyStateSpace(ref lamb, ref Vec, num);
                Normalize(Vec, MM, massNormalize, normalize);
                return lamb;
            }
            return Eig_(KK, MM, out Vec, sort, num, massNormalize, normalize, ascending, sortMode);
        }
        /// <summary>
        /// Added by MRB
        /// </summary>
        /// <param name="K">stiffeness matrix</param>
        /// <param name="Vec">eigen vectors, by columns</param>
        /// <param name="sort">sort by values</param>
        /// <param name="ascending">ascending sort?</param>
        /// <param name="num">number of eigen values and eigen vectors</param>
        /// <param name="normalize">vectors will be normalized</param>
        /// <param name="massNormalize">vectors will be normalized due to given mass</param>
        /// <param name="mass">will be used for vector normalizition and valuse will by divide by mass</param>
        /// <returns>eigen values of K divided by mass</returns>
        public static Vector Eig(Matrix K, out Matrix Vec, bool sort = true, bool ascending = true, int num = -1, bool normalize = false, bool massNormalize = false, double mass = 1.0)
        {
            if (num <= 0)
                num = K.ColumnCount;
            Vector Val = null;
            Matrix Vec0;
            if (ARMADILO.Enabled)
            {
                ComplexMatrix Vec_;
                var Val_ = ARMADILO.eig_gen(K, out Vec_);
                Val = Val_.Real();
                Vec0 = Vec_.Real();
            }
            else
            {
                var Val_ = (new EigenSystem()).GetEigenvalues_Real(K, out Vec0);
                Val = Val_.GetColumnVector(0);
            }
            if (sort)
            {
                for (int i = 0; i < num; i++)
                {
                    for (int j = i + 1; j < Val.Count; j++)
                    {
                        if (ascending)
                        {
                            if (Val[i] > Val[j])
                            {
                                Val.Swap(i, j);
                                Vec0.SwapColumns(i, j);
                            }
                        }
                        else
                            if (Val[i] < Val[j])
                        {
                            Val.Swap(i, j);
                            Vec0.SwapColumns(i, j);
                        }
                    }
                    if (normalize)
                    {
                        var V = Vec0.GetColumnVector(i);
                        V = V.Normalize();
                        Vec0.InsertAt(0, i, V, -1, 1);
                    }
                    else if (massNormalize)
                    {
                        var V = Vec0.GetColumnVector(i);
                        V = V.Normalize() * Math.Sqrt(mass);
                        Vec0.InsertAt(0, i, V, -1, 1);
                    }
                }
            }
            Vec = new Matrix(K.RowCount, num);
            Vec.InsertAt(0, 0, Vec0, colCount: num);
            Val.Resize(num);
            Val /= mass;
            return Val;
        }


        /// <summary>
        /// Added by MRB
        /// </summary>
        /// <param name="K">stiffeness matrix</param>
        /// <param name="Vec">eigen vectors, by columns</param>
        /// <param name="sort">sort by values</param>
        /// <param name="ascending">ascending sort?</param>
        /// <param name="num">number of eigen values and eigen vectors</param>
        /// <param name="normalize">vectors will be normalized</param>
        /// <param name="massNormalize">vectors will be normalized due to given mass</param>
        /// <param name="mass">will be used for vector normalizition and valuse will by divide by mass</param>
        /// <returns>eigen values of K divided by mass</returns>
        public static Vector Eig(SymmetricMatrix K, out Matrix Vec, bool sort = true, bool ascending = true, int num = -1, bool normalize = false, bool massNormalize = false, double mass = 1.0)
        {
            if (num <= 0)
                num = K.ColumnCount;

            Vector Val = null;
            Matrix Vec0;
            if (ARMADILO.Enabled)
            {
                ComplexMatrix Vec_;
                var Val_ = ARMADILO.eig_gen(K, out Vec_);
                Val = Val_.Real();
                Vec0 = Vec_.Real();
            }
            else
            {
                var Val_ = (new EigenSystem()).GetEigenvalues_Real(K, out Vec0);
                Val = Val_.GetColumnVector(0);
            }
            if (sort)
            {
                for (int i = 0; i < num; i++)
                {
                    for (int j = i + 1; j < Val.Count; j++)
                    {
                        if (ascending)
                        {
                            if (Val[i] > Val[j])
                            {
                                Val.Swap(i, j);
                                Vec0.SwapColumns(i, j);
                            }
                        }
                        else
                            if (Val[i] < Val[j])
                        {
                            Val.Swap(i, j);
                            Vec0.SwapColumns(i, j);
                        }
                    }
                    if (normalize)
                    {
                        var V = Vec0.GetColumnVector(i);
                        V = V.Normalize();
                        Vec0.InsertAt(0, i, V, -1, 1);
                    }
                    else if (massNormalize)
                    {
                        var V = Vec0.GetColumnVector(i);
                        V = V.Normalize() * Math.Sqrt(mass);
                        Vec0.InsertAt(0, i, V, -1, 1);
                    }
                }
            }
            Vec = new Matrix(K.RowCount, num);
            Vec.InsertAt(0, 0, Vec0, colCount: num);
            Val.Resize(num);
            Val /= mass;
            return Val;
        }


        /// <summary>
        /// A[n]*x^(n) + ... + A[2]*x_dot_dot + A[1]*x_dot + A[0]*x = 0,
        /// n = A.Count-1
        /// </summary>
        /// <param name="A"></param>
        /// <param name="EigenVec"></param>
        /// <returns></returns>
        public static ComplexVector Eig(List<ComplexMatrix> A, out ComplexMatrix EigenVec,
            int num = -1, bool massNormalize = true,
            bool normalize = false, bool sort = true)
        {
            if (!ARMADILO.Enabled)
                throw new NotImplementedException("ComplexMatrix.EigenValue: only with armadillo");
            var n = A.Count - 1;
            var rows = A[0].RowCount;
            var M = new ComplexMatrix(n * rows);
            var K = new ComplexMatrix(n * rows);
            K.InsertAt(0, 0, A[0]);
            for (int i = 0; i < n; i++)
                M.InsertAt(0, i * rows, -1 * A[i + 1]);
            var I = Matrix.I(rows);
            for (int i = 1; i < n; i++)
            {
                M.InsertAt(i * rows, (i - 1) * rows, I);
                K.InsertAt(i * rows, i * rows, I);
            }
            ComplexMatrix Vec;
            var vals = ARMADILO.eig_pair(K, M, out Vec);
            if (massNormalize || normalize)
                Normalize(Vec, M, massNormalize, normalize);
            if (sort)
            {
                //SortFrequencyStateSpace(ref vals, ref Vec, -1);

                for (int i = 0; i < vals.Count; i++)
                {
                    for (int j = i + 1; j < vals.Count; j++)
                        if (Vec.GetColumnVector(i).MaxAbs(0, 2 * rows) < Vec.GetColumnVector(j).MaxAbs(0, 2 * rows))
                        {
                            vals.Swap(i, j);
                            Vec.SwapColumns(i, j);
                        }
                }
            }
            else
                rows = -1;
            if (num > vals.Count) num = vals.Count;
            EigenVec = Vec[0, 0, rows, num];
            vals = vals[0, num];
            return vals;
        }

        public static void SortByTracking(
           List<double> x_axis,
           List<ComplexVector> values,
           List<ComplexMatrix> vectors = null,
           int reference = 0)
        {
            if (values.Count <= 1) return;
            var n = values[0].Count;
            for (int i = reference + 1; i < values.Count; i++)
            {
                var v_ref = values[reference];
                if (i - reference > 1)
                {
                    var v1 = values[i - 2];
                    var v2 = values[i - 1];
                    var x1 = x_axis[i - 2];
                    var x2 = x_axis[i - 1];
                    var x3 = x_axis[i];
                    v_ref = v2 + (v2 - v1) * ((x3 - x2) / (x2 - x1));
                }
                for (var j = 0; j < n; j++)
                {
                    var k = (values[i] - v_ref[j]).MinAbsIndex();
                    values[i].Swap(k, j);
                    if (vectors != null)
                        vectors[i].SwapColumns(k, j);
                }
            }
            for (int i = reference - 1; i >= 0; i--)
            {
                var v_ref = values[reference];
                if (reference - i > 1)
                {
                    var v1 = values[i + 2];
                    var v2 = values[i + 1];
                    var x1 = x_axis[i + 2];
                    var x2 = x_axis[i + 1];
                    var x3 = x_axis[i];
                    v_ref = v2 + (v2 - v1) * ((x3 - x2) / (x2 - x1));
                }
                for (var j = 0; j < n; j++)
                {
                    var k = (values[i] - v_ref[j]).MinAbsIndex();
                    values[i].Swap(k, j);
                    if (vectors != null)
                        vectors[i].SwapColumns(k, j);
                }
            }
        }


        #endregion

    }
}
