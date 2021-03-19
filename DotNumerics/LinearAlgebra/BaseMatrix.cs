#region Copyright � 2009, De Santiago-Castillo JA. All rights reserved.

//Copyright � 2009 Jose Antonio De Santiago-Castillo 
//E-mail:JAntonioDeSantiago@gmail.com
//Web: www.DotNumerics.com
//
#endregion

using System;
using System.Collections.Generic;
using System.Text;
using System.Diagnostics;
using System.IO;
using DotNumerics.LinearAlgebra.CSLapack;
using System.Threading;
using System.ComponentModel;
#if NET4
using System.Threading.Tasks;
#endif


namespace DotNumerics.LinearAlgebra
{

    /// <summary>
    /// Represents a Base Matrix.
    /// </summary>
    [DebuggerDisplay("[{RowCount},{ColumnCount}]")]
    //[DebuggerDisplay("[{RowCount},{ColumnCount}]", Name = "BandMatrix")]
    [DebuggerTypeProxy(typeof(MatrixDebuggerDisplay))]
    public abstract class BaseMatrix : IMatrix<double>
    {
        public static bool AllowEmptyMatrices = false;

        public delegate void OnChangedEvent(BaseMatrix sender, int row, int col);
        public OnChangedEvent onValueChanged = null;

        #region Static Fields

        internal static DGETRF _dgetrf;
        internal static DGETRI _dgetri;

        #endregion


        #region Fields
        /// <summary>
        /// Los datos de la matriz, los datos se almacenan en un un array unidimensional,
        /// Los elementos se almacenan por columnas, esto para que sean compatible con los Arrays de Fortran
        /// </summary>
        [DebuggerBrowsable(DebuggerBrowsableState.Never)]
        protected double[] _Data;
        /// <summary>
        /// El numero de renglones
        /// </summary>
        [DebuggerBrowsable(DebuggerBrowsableState.Never)]
        protected int _RowCount;

        /// <summary>
        /// El numero de columnas
        /// </summary>
        [DebuggerBrowsable(DebuggerBrowsableState.Never)]
        protected int _ColumnCount;

        #endregion


        #region  Public Constructors

        /// <summary>
        /// Initializes a new instance of the BaseMatrix class of the given size.
        /// </summary>
        /// <param name="rows">Number of rows.</param>
        /// <param name="columns">Number of columns.</param>
        public BaseMatrix(int rows, int columns)
        {
            if (!AllowEmptyMatrices)
            {
                if (rows < 1) throw new System.ArgumentException("rows < 1");
                if (columns < 1) throw new System.ArgumentException("columns < 1");
            }

            this._Data = new double[rows * columns];
            this._RowCount = rows;
            this._ColumnCount = columns;
        }

        /// <summary>
        /// Initializes a new instance of the BaseMatrix class of the given size using a array
        /// </summary>
        /// <param name="rows">Number of rows.</param>
        /// <param name="columns">Number of columns.</param>
        /// <param name="data">The data, the data is copied.</param>
        internal BaseMatrix(int rows, int columns, double[] data, bool MakeACopyOfData = true)
        {
            if (!AllowEmptyMatrices)
            {
                if (rows < 1) throw new System.ArgumentException("rows < 1");
                if (columns < 1) throw new System.ArgumentException("columns < 1");
            }
            this._RowCount = rows;
            this._ColumnCount = columns;
            if (MakeACopyOfData)
            {
                this._Data = new double[rows * columns];
                data.CopyTo(this._Data, 0);
            }
            else
                this._Data = data;
            ////Si incluye la posibilidad de que los datos tengan menos valores que la matriz a crear 
            //for (int i = 0; i < Math.Min(this.MeData.Length, data.Length); i++)
            //{
            //    this.MeData[i] = data[i];
            //}
        }

        /// <summary>
        /// Initializes a new instance of the BaseMatrix class of the given size.
        /// </summary>
        /// <param name="size">Size</param>
        public BaseMatrix(int size)
        {
            if (!AllowEmptyMatrices)
                if (size < 1) throw new System.ArgumentException("size < 1");

            this._Data = new double[size * size];
            this._RowCount = size;
            this._ColumnCount = size;
        }

        /// <summary>
        /// Initializes a new instance of the BaseMatrix class of the given size using a array
        /// </summary>
        /// <param name="size">Size</param>
        /// <param name="data">The data</param>
        internal BaseMatrix(int size, double[] data)
        {
            if (!AllowEmptyMatrices)
                if (size < 1) throw new System.ArgumentException("size < 1");

            this._Data = new double[size * size];
            this._RowCount = size;
            this._ColumnCount = size;

            data.CopyTo(this._Data, 0);

            ////Si incluye la posibilidad de que los datos tengan menos valores que la matriz a crear 
            //for (int i = 0; i < Math.Min(this.MeData.Length, data.Length); i++)
            //{
            //    this.MeData[i] = Data[i];
            //}
        }

        internal BaseMatrix(int rows, int columns, Func<int, int, double> filler) : this(rows, columns)
        {
            Fill(filler);
        }
        internal BaseMatrix(int rows, int columns, Func<int, double> filler) : this(rows, columns)
        {
            Fill_Diag(filler);
        }

        public void Fill(Func<int, int, double> filler)
        {
            Parallel.For(0, ColumnCount, j =>
            {
                var jj = j * this._RowCount;
                for (int i = 0; i < RowCount; i++)
                    Data[jj + i] = filler(i, j);
            });
        }
        public void Fill_Diag(Func<int, double> filler)
        {
            for (int i = 0; i < Math.Min(RowCount, ColumnCount); i++)
                this[i, i] = filler(i);
        }
        #endregion


        #region Public Properties

        /// <summary>
        /// Los datos de la matriz
        /// </summary>
        [DebuggerBrowsable(DebuggerBrowsableState.Never)]
        public double[] Data
        {
            get { return this._Data; }
        }
        public void SetData(double[] Data)
        {
            _Data = Data;
        }

        /// <summary>
        /// Returns the number of rows.
        /// </summary>
        [DebuggerBrowsable(DebuggerBrowsableState.Never)]
        public int RowCount
        {
            get { return _RowCount; }
            //set { MeRowCount = value; }
        }
        /// <summary>
        /// Returns the number of columns.
        /// </summary>
        [DebuggerBrowsable(DebuggerBrowsableState.Never)]
        public int ColumnCount
        {
            get { return _ColumnCount; }
            //set { MeColumnCount = value; }
        }

        /// <summary>
        /// Gets a value indicating if the matrix is square.
        /// </summary>
        public bool IsSquare
        {
            get
            {
                bool isSquare = false;
                if (this._ColumnCount == this.RowCount) isSquare = true;
                return isSquare;
            }
        }
        /// <summary>
        /// Gets or set the value of a element of this matrix.
        /// </summary>
        /// <param name="row">The row value (zero-based).</param>
        /// <param name="column">The column value (zero-based).</param>
        /// <returns>The matrix element at (row, column).</returns>
        public virtual double this[int row, int column]
        {
            get
            {
                if (column >= this._ColumnCount)
                {
                    throw new ArgumentException("Index was outside the bounds of the matrix (get). [" + row + "," + column + "]>[" + RowCount + "," + ColumnCount + "]");
                }

                return this._Data[row + column * this._RowCount];
            }
            set
            {
                if (column >= this._ColumnCount)
                {
                    throw new ArgumentException("Index was outside the bounds of the matrix (set). [" + row + "," + column + "]>[" + RowCount + "," + ColumnCount + "]");
                }
                if (onValueChanged != null)
                    onValueChanged(this, row, column);

                this._Data[row + column * this._RowCount] = value;
            }
        }
        /// <summary>
        /// sub matrix
        /// </summary>
        /// <param name="startRow"></param>
        /// <param name="startCol"></param>
        /// <param name="rowCount"></param>
        /// <param name="colCount"></param>
        /// <returns></returns>
        public Matrix this[int startRow, int startCol, int rowCount, int colCount]
        {
            get { return SubMatrix(startRow, startCol, rowCount, colCount); }
            set { InsertAt(startRow, startCol, value, rowCount, colCount); }
        }
        public Matrix this[List<int> Rows, List<int> Columns]
        {
            get { return SubMatrix(Rows, Columns); }
            set { InsertAt(Rows, Columns, value); }
        }

        public double xx
        {
            get { return this[0, 0]; }
            set { this[0, 0] = value; }
        }
        public double xy
        {
            get { return this[0, 1]; }
            set { this[0, 1] = value; }
        }
        public double xz
        {
            get { return this[0, 2]; }
            set { this[0, 2] = value; }
        }
        public double yx
        {
            get { return this[1, 0]; }
            set { this[1, 0] = value; }
        }
        public double yy
        {
            get { return this[1, 1]; }
            set { this[1, 1] = value; }
        }
        public double yz
        {
            get { return this[1, 2]; }
            set { this[1, 2] = value; }
        }
        public double zx
        {
            get { return this[2, 0]; }
            set { this[2, 0] = value; }
        }
        public double zy
        {
            get { return this[2, 1]; }
            set { this[2, 1] = value; }
        }
        public double zz
        {
            get { return this[2, 2]; }
            set { this[2, 2] = value; }
        }

        #endregion

        #region	 Private Methods

        /// <summary>Check if size(this) == size(B) </summary>
        internal protected void CheckMatrixDimensions(BaseMatrix B)
        {
            if (this._RowCount != B.RowCount || B.ColumnCount != this._ColumnCount)
            {
                throw new System.ArgumentException("Matrix dimensions must agree.");
            }
        }
        #endregion //  Private Methods

        #region Elementary linear operations

        /// <summary>
        /// aij=Math.Abs(aij)
        /// </summary>
        public virtual void ElementsAbs()
        {
            for (int i = 0; i < this._Data.Length; i++)
            {
                this._Data[i] = Math.Abs(this._Data[i]);
            }
        }

        /// <summary>
        /// Element-by-element division: aij = aij/bij
        /// </summary>
        /// <param name="B">The B Matrix.</param>
        public virtual void ElemntsDiv(BaseMatrix B)
        {
            CheckMatrixDimensions(B);
            double[] BData = B.Data;
            for (int i = 0; i < this._Data.Length; i++)
            {
                this._Data[i] /= BData[i];
            }
        }

        /// <summary>
        /// Element-by-element multiplication: aij = aij*bij
        /// </summary>
        /// <param name="B">The B Matrix.</param>
        public virtual void ElemntsMult(BaseMatrix B)
        {
            CheckMatrixDimensions(B);
            double[] BData = B.Data;
            Parallel.For(0, this._Data.Length, i =>
            {
                this._Data[i] *= BData[i];
            });
        }

        /// <summary>
        /// Addition C=A+B
        /// </summary>
        /// <param name="B">The Matrix</param>
        /// <returns>C=A+B</returns>
        public virtual Matrix Add(BaseMatrix B, int thread_count = -1)
        {

            CheckMatrixDimensions(B);

            Matrix C = new Matrix(this._RowCount, this._ColumnCount);

            double[] BData = B.Data;
            double[] dataC = C.Data;
            Parallel.For(0, _Data.Length, i =>
            {
                dataC[i] = this._Data[i] + BData[i];
            });
            return C;
        }




        /// <summary>
        /// In place scalar-matrix multiplication, A=s*A
        /// </summary>
        /// <param name="s">The scalar</param>
        public virtual void MultiplyInplace(double s)
        {
            Parallel.For(0, _Data.Length, i =>
            {
                this._Data[i] *= s;
            });
        }

        /// <summary>
        /// Scalar-matrix multiplication, C=s*A
        /// </summary>
        /// <param name="s">The scalar</param>
        /// <returns>C=s*A</returns>
        public Matrix Multiply(double s, int thread_count = -1)
        {
            Matrix C = new Matrix(this._RowCount, this._ColumnCount);
            double[] dataC = C.Data;
            Parallel.For(0, _Data.Length, i =>
            {
                dataC[i] = this._Data[i] * s;
            });
            return C;
        }
        public ComplexMatrix Multiply(Complex s, int thread_count = -1)
        {
            var C = new ComplexMatrix(this._RowCount, this._ColumnCount);
            var dataC = C.Data;
            Parallel.For(0, _Data.Length, i =>
            {
                dataC[i] = this._Data[i] * s;
            });
            return C;
        }

        /// <summary>
        /// Matrix-Matrix multiplication, C=A*B
        /// </summary>
        /// <param name="B">The matrix.</param>
        /// <returns>C=A*B</returns>
        public Matrix Multiply(BaseMatrix B, int thread_count = -1)
        {
            if (B.RowCount != this.ColumnCount)
            {
                throw new System.ArgumentException("Matrix dimensions are not valid.");
            }

            if (ARMADILO.Enabled && this is Matrix && B is Matrix)
                return ARMADILO.mat_multiply((Matrix)this, (Matrix)B);

            Matrix C = new Matrix(this.RowCount, B.ColumnCount);

            double[] AData = this.Data;
            double[] BData = B.Data;
            double[] CData = C.Data;

            int ARows = this.RowCount;
            int AColumns = this.ColumnCount;

            int BRows = B.RowCount;
            int BColumns = B.ColumnCount;

            Parallel.For(0, BColumns, j =>
                 {
                     int indexBJ = j * BRows;
                     int indexAJ = j * ARows;
                     for (int i = 0; i < ARows; i++)
                     {
                         var Sum = 0.0;
                         for (int k = 0; k < AColumns; k++)
                         {
                             Sum += AData[i + k * ARows] * BData[k + indexBJ];
                         }
                         CData[i + indexAJ] = Sum;
                     }
                 });
            return C;
            #region ?
            //To reading every time elements from array , why we are taking some group of element i.e. Block size, then no need to read every element. A groups of element will be on catche and we can do fast as given above algo. This algorithm called " Block Algorithm". This Block algorithm can be applied many place where this type of situation will come.

            //Block Algorithm for Matrix Multiplication:

            //Code: C
            //            #define n 1000
            //int main()
            //{
            //    int a[n][n],b[n][n],c[n][n];
            //    c[0][0]=0;
            //    for( i=0;i<n;++i)
            //    {
            //        for(j=0;j<n;++j)
            //        {
            //            for(k=0;k<n;++k)
            //            {
            //                c[i][j] = c[i][j] + a[i][k] * b[k][j]
            //            }
            //        }
            //    }
            //    return 0;
            //}
            //To reading every time elements from array , why we are taking some group of element i.e. Block size, then no need to read every element. A groups of element will be on catche and we can do fast as given above algo. This algorithm called " Block Algorithm". This Block algorithm can be applied many place where this type of situation will come.

            //Block Algorithm for Matrix Multiplication:

            //Code: C

            //#define n 1000
            //#define BlockSize  100
            //int main()
            //{
            //    int a[n][n],b[n][n],c[n][n];
            //    c[0][0]=0;
            //    for( i1=0;i1<(n/BlockSize);++i1)
            //    {
            //        for(j1=0;j1<(n/BlockSize);++j1)
            //        {
            //            for(k1=0;k1<(n/BlockSize);++k1)
            //            {
            //                for(i=i1=0;i<min(i1+BlockSize-1);++i)
            //                {
            //                    for(j=j1=0;j<min(j1+BlockSize-1);++j)
            //                    {
            //                        for(k=k1;k<min(k1+BlockSize-1);++k)
            //                        {               
            //                            c[i][j] = c[i][j] + a[i][k] * b[k][j]
            //                        }
            //                    }
            //                }
            //            }
            //        }
            //           }
            // return 0;
            //}
            #endregion
        }

        /// <summary>
        /// Matrix subtraction, C=A-B
        /// </summary>
        /// <param name="B">The Matrix</param>
        /// <returns>C=A-B</returns>
        public Matrix Subtract(BaseMatrix B, int thread_count = -1)
        {

            CheckMatrixDimensions(B);

            Matrix C = new Matrix(this._RowCount, this._ColumnCount);

            double[] BData = B.Data;
            double[] dataC = C.Data;
            Parallel.For(0, this._Data.Length, i =>
               {
                   dataC[i] = this._Data[i] - BData[i];
               });
            return C;
        }

        /// <summary>
        /// In place unary minus -A
        /// </summary>
        public virtual void UnaryMinusInplace(int thread_count = -1)
        {
            Parallel.For(0, this._Data.Length, i =>
            {
                this._Data[i] = -this._Data[i];
            });
        }

        #endregion

        #region Methods

        public string ToString1(string format = "", bool simple = true, int digits = 4)
        {
            var res = "";
            if (!simple)
                res += RowCount + "*" + ColumnCount + "\r\n----------\r\n";
            int max = 0;
            for (int r = 0; r < RowCount; r++)
            {
                for (int c = 0; c < ColumnCount; c++)
                {
                    if (format == "")
                    {
                        var s = this[r, c].ToString("G" + digits);
                        max = Math.Max(max, s.Length);
                        res += '>' + s.PadLeft(digits + 6, ' ');
                    }
                    else
                        res += (this[r, c] >= 0 ? " " : "") + this[r, c].ToString(format);
                    if (!simple)
                        res += ", ";
                    else
                        res += "  ";
                }
                if (!simple)
                    res += "  ;\r\n";
                else
                    res += "\r\n";
            }
            if (format == "")
            {
                var s = ">".PadRight(digits + 6 - max + 1, ' ');
                res = res.Replace(s, "");
            }
            return res;
        }

        /// <summary>
        /// Resizes the matrix
        /// Added by MRB
        /// </summary>
        /// <param name="newRowCount"></param>
        /// <param name="newColumnCount"></param>
        public void Resize(int newRowCount, int newColumnCount, bool keep_elems = true)
        {
            if (keep_elems && newRowCount != this.RowCount)
            {
                var tmp = new Matrix(this._RowCount, this._ColumnCount);
                Array.Copy(_Data, tmp._Data, _Data.Length);
                _Data = new double[newRowCount * newColumnCount];
                this._RowCount = newRowCount;
                this._ColumnCount = newColumnCount;

                Parallel.For(0, Math.Min(tmp.RowCount, this.RowCount), i =>
                {
                    for (int j = 0; j < Math.Min(tmp.ColumnCount, this.ColumnCount); j++)
                        this[i, j] = tmp[i, j];
                });
            }
            else
            {
                Array.Resize(ref _Data, newRowCount * newColumnCount);
                this._RowCount = newRowCount;
                this._ColumnCount = newColumnCount;
            }
        }

        /// <summary>
        /// Calculates the inverse of the matrix.
        /// </summary>
        /// <returns>The inverse of the matrix.</returns>
        public Matrix Inverse()
        {

            if (this.IsSquare != true)
                throw new System.ArgumentException("This is not a square matrix. (hint: use PseudoInverse)");

            if (ARMADILO.Enabled && this is Matrix)
                return ARMADILO.inv((Matrix)this);

            if (BaseMatrix._dgetrf == null)
            {
                BaseMatrix._dgetrf = new DGETRF();
            }
            if (BaseMatrix._dgetri == null)
            {
                BaseMatrix._dgetri = new DGETRI();
            }


            Matrix inverseMatrix = new Matrix(this.RowCount, this.ColumnCount, this.Data);

            double[] inverseData = inverseMatrix.Data;

            int[] ipiv = new int[this.RowCount];



            int Info = 0;

            double[] Work = new double[1];
            int LWork = -1;

            //Calculamos LWORK 
            BaseMatrix._dgetri.Run(this.RowCount, ref inverseData, 0, this.RowCount, ipiv, 0, ref Work, 0, LWork, ref Info);
            LWork = Convert.ToInt32(Work[0]);
            if (LWork > 0)
            {
                Work = new double[LWork];

                BaseMatrix._dgetrf.Run(this.RowCount, this.ColumnCount, ref inverseData, 0, this.RowCount, ref ipiv, 0, ref Info);


                #region Error
                /// = 0:  successful exit
                /// .LT. 0:  if INFO = -i, the i-th argument had an illegal value
                /// .GT. 0:  if INFO = i, U(i,i) is exactly zero. The factorization
                /// has been completed, but the factor U is exactly
                /// singular, and division by zero will occur if it is used
                /// to solve a system of equations.

                if (Info < 0)
                {
                    string infoSTg = Math.Abs(Info).ToString();
                    throw new ArgumentException("the " + infoSTg + " -th argument had an illegal value");
                }
                else if (Info > 0)
                {
                    string infoSTg = Math.Abs(Info).ToString();
                    throw new Exception("The matrix is numerically singular..");
                }

                #endregion


                BaseMatrix._dgetri.Run(this.RowCount, ref inverseData, 0, this.RowCount, ipiv, 0, ref Work, 0, LWork, ref Info);
            }
            else
            {

                //Error
            }

            #region Error
            /// (output) INTEGER
            /// = 0:  successful exit
            /// .LT. 0:  if INFO = -i, the i-th argument had an illegal value
            /// .GT. 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
            /// singular and its inverse could not be computed.

            if (Info < 0)
            {
                string infoSTg = Math.Abs(Info).ToString();
                throw new ArgumentException("the " + infoSTg + " -th argument had an illegal value");
            }
            else if (Info > 0)
            {
                string infoSTg = Math.Abs(Info).ToString();
                throw new Exception("The matrix is numerically singular..");
            }

            #endregion

            return inverseMatrix;
        }


        /// <summary>
        /// Calculates the determinant of the matrix.
        /// </summary>
        /// <returns>The determinant of the matrix.</returns>
        public double Determinant()
        {
            double det = 1.0;


            if (this.IsSquare != true)
            {
                throw new System.ArgumentException("This is not a square matrix.");
            }
            if (ARMADILO.Enabled && this is Matrix)
                return ARMADILO.det(this as Matrix);

            if (BaseMatrix._dgetrf == null)
            {
                BaseMatrix._dgetrf = new DGETRF();
            }


            Matrix clonMatrix = new Matrix(this.RowCount, this.ColumnCount, this.Data);

            double[] clonData = clonMatrix.Data;

            int[] ipiv = new int[this.RowCount];



            int Info = 0;



            BaseMatrix._dgetrf.Run(this.RowCount, this.ColumnCount, ref clonData, 0, this.RowCount, ref ipiv, 0, ref Info);


            #region Error
            // <param name="INFO">
            // (output) INTEGER
            //= 0:  successful exit
            // .LT. 0:  if INFO = -i, the i-th argument had an illegal value
            // .GT. 0:  if INFO = i, U(i,i) is exactly zero. The factorization
            // has been completed, but the factor U is exactly
            // singular, and division by zero will occur if it is used
            // to solve a system of equations.
            if (Info < 0)
            {
                string infoSTg = Math.Abs(Info).ToString();
                throw new ArgumentException("the " + infoSTg + " -th argument had an illegal value");
            }
            //else if (Info > 0)
            //{
            //    string infoSTg = Math.Abs(Info).ToString();
            //    throw new Exception("The matrix is numerically singular..");
            //}

            #endregion



            //The LU factorization yields three matrices, the product of which is the original 
            //complex matrix. Therefore the determinat is the product of the three determinants 
            //of P, L and U. The determinant of the triangular matrices L and U is the product 
            //of the elements on the diagonal - as for any triangular matrix (for L this is 1 
            //as all elements of the diagonal are one.) The determinant of P is either +1 or -1 
            //depending of whether the number of row permutations is even or odd. 

            //Thank you very much for your answer. It seems to be that your message got tructated somehow, but I think I got the point.
            //I also did some searching on the web and found the following two pieces of code which both claim to calculate the determinant of a square matrix.

            //============================================

            //call dgetrf(n,n,a,n,piv,info)
            //det = 0d0
            //if (info.ne.0) then
            //return
            //endif
            //det = 1d0
            //do 10,i=1,n
            //if (piv(i).ne.i) then
            //det = -det * a(i,i)
            //else
            //det = det * a(i,i)
            //endif
            //10 continue
            //end



            for (int i = 0; i < this._RowCount; i++)
            {
                if (ipiv[i] != (i + 1))  // i+1 debido a que aqui la base es 0 y en fortran es 1
                {
                    det *= -clonMatrix[i, i];
                }
                else
                {
                    det *= clonMatrix[i, i];
                }

            }

            return det;
        }



        /// <summary>
        ///  Gets a column vector of this matrix at the selected position.
        /// </summary>
        /// <param name="columnIndex">The column index (zero-based).</param>
        /// <returns>The column vector.</returns>
        public Vector GetColumnVector(int columnIndex)
        {

            if (columnIndex >= this._ColumnCount)
            {
                throw new System.ArgumentException("columnIndex >= number of columns.");
            }

            if (columnIndex < 0)
            {
                throw new System.ArgumentException("columnIndex < 0");
            }

            var columnVect = new Vector(VectorType.Column, this._RowCount);
            Array.Copy(this._Data, columnIndex * this._RowCount, columnVect.Data, 0, this._RowCount);
            return columnVect;
        }

        /// <summary>
        /// Gets the column vectors of this matrix.
        /// </summary>
        /// <returns>The columns vectors.</returns>
        public Vector[] GetColumnVectors()
        {
            var res = new Vector[ColumnCount];
            for (int i = 0; i < ColumnCount; i++)
            {
                var columnVect = new Vector(VectorType.Column, this._RowCount);
                Array.Copy(this._Data, i * this._RowCount, columnVect.Data, 0, this._RowCount);
                res[i] = columnVect;
            }
            return res;
        }

        /// <summary>
        ///  Gets a column array of this matrix at the selected position.
        /// </summary>
        /// <param name="columnIndex">The column index (zero-based).</param>
        /// <returns>The column array.</returns>
        public double[] GetColumnArray(int columnIndex)
        {

            if (columnIndex >= this._ColumnCount)
            {
                throw new System.ArgumentException("columnIndex >= number of columns.");
            }

            if (columnIndex < 0)
            {
                throw new System.ArgumentException("columnIndex < 0");
            }

            double[] VectData = new double[this._RowCount];
            for (int i = 0; i < VectData.Length; i++)
            {
                VectData[i] = this._Data[i + columnIndex * this._RowCount];
            }

            return VectData;
        }


        /// <summary>
        /// Gets the row vectors of this matrix.
        /// </summary>
        /// <returns>The row vectors.</returns>
        public Vector[] GetRowVectors()
        {
            Vector[] rowVects = new Vector[this.RowCount];

            double[] VectData;
            for (int i = 0; i < this._RowCount; i++)
            {
                rowVects[i] = new Vector(VectorType.Row, this._ColumnCount);
                VectData = rowVects[i].Data;
                for (int j = 0; j < VectData.Length; j++)
                {
                    VectData[j] = this._Data[i + j * this._RowCount];
                }
            }

            return rowVects;
        }


        /// <summary>
        ///  Gets a row vector of this matrix at the selected position.
        /// </summary>
        /// <param name="rowIndex">The row index (zero-based).</param>
        /// <returns>The row vector.</returns>
        public Vector GetRowVector(int rowIndex)
        {

            if (rowIndex >= this._RowCount)
            {
                throw new System.ArgumentException("rowIndex >= number of rows.");
            }

            if (rowIndex < 0)
            {
                throw new System.ArgumentException("rowIndex < 0");
            }

            Vector rowVect;

            double[] VectData;
            rowVect = new Vector(VectorType.Row, this._ColumnCount);
            VectData = rowVect.Data;
            for (int j = 0; j < VectData.Length; j++)
            {
                VectData[j] = this._Data[rowIndex + j * this._RowCount];
            }

            return rowVect;
        }

        /// <summary>
        ///  Gets a row array of this matrix at the selected position.
        /// </summary>
        /// <param name="rowIndex">The row index (zero-based).</param>
        /// <returns>The row array.</returns>
        public double[] GetRowArray(int rowIndex)
        {

            if (rowIndex >= this._RowCount)
            {
                throw new System.ArgumentException("rowIndex >= number of rows.");
            }

            if (rowIndex < 0)
            {
                throw new System.ArgumentException("rowIndex < 0");
            }

            double[] VectData = new double[this._ColumnCount];

            for (int j = 0; j < VectData.Length; j++)
            {
                VectData[j] = this._Data[rowIndex + j * this._RowCount];
            }
            //}

            return VectData;
        }


        /// <summary>
        /// Returns the equivalent string representation of the matrix.
        /// </summary>
        /// <returns>The string representation of the  matrix.</returns>
        public string MatrixToString()
        {
            using (StringWriter writer = new StringWriter())
            {
                for (int i = 0; i < this._RowCount; i++)
                {
                    for (int j = 0; j < this._ColumnCount; j++)
                        writer.Write(this[i, j].ToString() + ", ");
                    writer.WriteLine();
                }
                return writer.ToString();
            }
        }

        /// <summary>
        /// Returns the equivalent string representation of the matrix.
        /// </summary>
        /// <param name="format">A numeric format string.</param>
        /// <returns>The string representation of the  matrix.</returns>
        public string MatrixToString(string format)
        {
            using (StringWriter writer = new StringWriter())
            {
                for (int i = 0; i < this._RowCount; i++)
                {
                    for (int j = 0; j < this._ColumnCount; j++)
                        writer.Write(this[i, j].ToString(format) + ", ");
                    writer.WriteLine();
                }
                return writer.ToString();
            }
        }

        /// <summary>
        /// Added by MRB
        /// </summary>
        /// <param name="i"></param>
        /// <param name="j"></param>
        public void SwapRows(int i, int j)
        {
            for (int c = 0; c < ColumnCount; c++)
            {
                var tmp = this[i, c];
                this[i, c] = this[j, c];
                this[j, c] = tmp;
            }
        }

        /// <summary>
        /// Added by MRB
        /// </summary>
        /// <param name="i"></param>
        /// <param name="j"></param>
        public void SwapColumns(int i, int j)
        {
            var i1 = _RowCount * i;
            var j1 = _RowCount * j;
            for (int c = 0; c < _RowCount; c++)
            {
                var tmp = _Data[i1 + c];
                _Data[i1 + c] = _Data[j1 + c];
                _Data[j1 + c] = tmp;
            }
            /*for (int r = 0; r < RowCount; r++)
            {
                var tmp = this[r, i];
                this[r, i] = this[r, j];
                this[r, j] = tmp;
            }*/
        }

        public void Split(int RowCount, int ColumnCount, ref Matrix A11, ref Matrix A12, ref Matrix A21, ref Matrix A22)
        {
            A11 = this[0, 0, RowCount, ColumnCount];
            A12 = this[0, ColumnCount, RowCount, -1];
            A21 = this[RowCount, 0, -1, ColumnCount];
            A22 = this[RowCount, ColumnCount, -1, -1];
        }
        public static Matrix Combine(Matrix A11, Matrix A12, Matrix A21, Matrix A22)
        {
            var A = new Matrix(A11.RowCount + A22.RowCount, A11.ColumnCount + A22.ColumnCount);
            A.InsertAt(0, 0, A11);
            A.InsertAt(0, A11.ColumnCount, A12);
            A.InsertAt(A11.RowCount, 0, A21);
            A.InsertAt(A11.RowCount, A11.ColumnCount, A22);
            return A;
        }
        public static Matrix CombineColumnMode(Matrix A1, Matrix A2)
        {
            var A = new Matrix(A1.RowCount + A2.RowCount, A1.ColumnCount);
            A.InsertAt(0, 0, A1);
            A.InsertAt(A1.RowCount, 0, A2);
            return A;
        }
        public static Matrix CombineRowMode(Matrix A1, Matrix A2)
        {
            var A = new Matrix(A1.RowCount, A1.ColumnCount + A2.ColumnCount);
            A.InsertAt(0, 0, A1);
            A.InsertAt(0, A1.ColumnCount, A2);
            return A;
        }
        public void Split(List<int> Rows1, List<int> Columns1, ref Matrix A11, ref Matrix A12, ref Matrix A21, ref Matrix A22)
        {
            var Rows2 = new List<int>();
            for (int i = 0; i < RowCount; i++)
                if (!Rows1.Contains(i))
                    Rows2.Add(i);
            var Columns2 = new List<int>();
            for (int i = 0; i < ColumnCount; i++)
                if (!Columns1.Contains(i))
                    Columns2.Add(i);
            A11 = A11 ?? new Matrix(Rows1.Count, Columns1.Count);
            A12 = A12 ?? new Matrix(Rows1.Count, Columns2.Count);
            A21 = A21 ?? new Matrix(Rows2.Count, Columns1.Count);
            A22 = A22 ?? new Matrix(Rows2.Count, Columns2.Count);

            for (int i = 0; i < Rows1.Count; i++)
            {
                for (int j = 0; j < Columns1.Count; j++)
                    A11[i, j] = this[Rows1[i], Columns1[j]];

                for (int j = 0; j < Columns2.Count; j++)
                    A12[i, j] = this[Rows1[i], Columns2[j]];
            }

            for (int i = 0; i < Rows2.Count; i++)
            {
                for (int j = 0; j < Columns1.Count; j++)
                    A21[i, j] = this[Rows2[i], Columns1[j]];

                for (int j = 0; j < Columns2.Count; j++)
                    A22[i, j] = this[Rows2[i], Columns2[j]];
            }
        }
        public static Matrix Cobine(List<int> Rows1, List<int> Columns1, Matrix A11, Matrix A12, Matrix A21, Matrix A22)
        {
            var A = new Matrix(A11.RowCount + A22.RowCount, A11.ColumnCount + A22.ColumnCount);
            var Rows2 = new List<int>();
            for (int i = 0; i < A.RowCount; i++)
                if (!Rows1.Contains(i))
                    Rows2.Add(i);
            var Columns2 = new List<int>();
            for (int i = 0; i < A.ColumnCount; i++)
                if (!Columns1.Contains(i))
                    Columns2.Add(i);

            for (int i = 0; i < Rows1.Count; i++)
            {
                for (int j = 0; j < Columns1.Count; j++)
                    A[Rows1[i], Columns1[j]] = A11[i, j];

                for (int j = 0; j < Columns2.Count; j++)
                    A[Rows1[i], Columns2[j]] = A12[i, j];
            }

            for (int i = 0; i < Rows2.Count; i++)
            {
                for (int j = 0; j < Columns1.Count; j++)
                    A[Rows2[i], Columns1[j]] = A21[i, j];

                for (int j = 0; j < Columns2.Count; j++)
                    A[Rows2[i], Columns2[j]] = A22[i, j];
            }
            return A;
        }



        /// <summary>
        /// One Norm for the matrix.
        /// </summary>
        /// <returns>The maximum column sum.</returns>
        public double Norm1()
        {
            double n = 0.0;
            double ColSum = 0.0;
            int NRows = this._RowCount;

            for (int j = 0; j < this._ColumnCount; j++)
            {
                ColSum = 0.0;
                for (int i = 0; i < this._RowCount; i++)
                {
                    ColSum += Math.Abs(this._Data[i + j * NRows]);

                }
                n = Math.Max(n, ColSum);
            }
            return n;
        }
        /// <summary>
        /// Infinity Norm for the matrix.
        /// </summary>
        /// <returns>The maximum row sum.</returns>
        public double NormInf()
        {
            double n = 0.0;
            double RowSum = 0.0;
            int NRows = this._RowCount;
            for (int i = 0; i < this._RowCount; i++)
            {
                RowSum = 0.0;
                for (int j = 0; j < this._ColumnCount; j++)
                {
                    RowSum += Math.Abs(this._Data[i + j * NRows]);
                }
                n = Math.Max(n, RowSum);
            }
            return n;
        }

        /// <summary>Frobenius norm</summary>
        /// <returns>The square root of sum of squares of all elements.</returns>
        public double FrobeniusNorm()
        {
            double n = 0;
            for (int i = 0; i < this._Data.Length; i++)
            {
                n = this.Hypot(n, this._Data[i]);
            }
            return n;
        }

        /// <summary>sqrt(a^2 + b^2) without under/overflow.</summary>
        private double Hypot(double a, double b)
        {
            double r;
            if (Math.Abs(a) > Math.Abs(b))
            {
                r = b / a;
                r = Math.Abs(a) * Math.Sqrt(1 + r * r);
            }
            else if (b != 0)
            {
                r = a / b;
                r = Math.Abs(b) * Math.Sqrt(1 + r * r);
            }
            else
            {
                r = 0.0;
            }
            return r;
        }

        /// <summary>
        /// Sum of elements =SUMij(A[i,j])
        /// </summary>
        /// <returns>The sum of elements.</returns>
        public double ElementsSum()
        {
            double TemSum = 0.0;
            for (int i = 0; i < this._Data.Length; i++)
            {
                TemSum += this._Data[i];
            }
            return TemSum;
        }
        /// <summary>
        /// Transposed matrix.
        /// </summary>
        /// <returns>The transposed matrix.</returns>
        public Matrix Transpose()
        {
            Matrix AT = new Matrix(this._ColumnCount, this._RowCount);
            int ATRows = AT.RowCount;
            int ATColumns = AT.ColumnCount;
            double[] ATData = AT.Data;
#if NET4
            Parallel.For(0, this._ColumnCount, j =>
            {
                for (int i = 0; i < this._RowCount; i++)
                    ATData[j + i * ATRows] = this._Data[i + j * this._RowCount];
            });
#else
            for (int j = 0; j < this._ColumnCount; j++)
                for (int i = 0; i < this._RowCount; i++)
                    ATData[j + i * ATRows] = this._Data[i + j * this._RowCount];
#endif

            return AT;
        }
        /// <summary>
        /// Transposed matrix.
        /// </summary>
        /// <returns>The transposed matrix.</returns>
        public Matrix T()
        {
            return Transpose();
        }

        /// <summary>Returns the trace of the matrix.</summary>
        /// <returns>Sum of the diagonal elements.</returns>
        public double Trace
        {
            get
            {
                double trace = 0;
                for (int i = 0; i < Math.Min(this.RowCount, this.ColumnCount); i++)
                {
                    trace += this[i, i];
                }
                return trace;
            }
        }


        /// <summary>
        /// Added By MRB
        /// </summary>
        /// <param name="startRow"></param>
        /// <param name="startCol"></param>
        /// <param name="M"></param>
        /// <param name="rowCount"></param>
        /// <param name="colCount"></param>
        public void InsertAt(int startRow, int startCol, Matrix M, int rowCount = -1, int colCount = -1)
        {
            if (rowCount == -1)
                rowCount = M.RowCount;
            if (colCount == -1)
                colCount = M.ColumnCount;
            Parallel.For(0, rowCount, i =>
            {
                for (int j = 0; j < colCount; j++)
                    this[i + startRow, j + startCol] = M[i, j];
            });
        }

        public void InsertAt(List<int> Rows, List<int> Columns, Matrix A1)
        {
            for (int i = 0; i < Rows.Count; i++)
                for (int j = 0; j < Columns.Count; j++)
                    this[Rows[i], Columns[j]] = A1[i, j];
        }

        /// <summary>
        /// Added By MRB
        /// </summary>
        /// <param name="startRow"></param>
        /// <param name="startCol"></param>
        /// <param name="M"></param>
        /// <param name="rowCount"></param>
        /// <param name="colCount"></param>
        public void PlusAt(int startRow, int startCol, Matrix M, int rowCount = -1, int colCount = -1)
        {
            if (rowCount == -1)
                rowCount = M.RowCount;
            if (colCount == -1)
                colCount = M.ColumnCount;
            Parallel.For(0, rowCount, i =>
            {
                for (int j = 0; j < colCount; j++)
                    this[i + startRow, j + startCol] += M[i, j];
            });
        }

        public Matrix SubMatrix(int startRow = 0, int startCol = 0, int rowCount = -1, int colCount = -1)
        {
            if (rowCount == -1)
                rowCount = RowCount - startRow;
            if (colCount == -1)
                colCount = ColumnCount - startCol;
            var res = new Matrix(rowCount, colCount);
            Parallel.For(0, rowCount, i =>
            {
                for (int j = 0; j < colCount; j++)
                    res[i, j] = this[i + startRow, j + startCol];
            });
            return res;
        }
        public Matrix SubMatrix(List<int> rows, List<int> columns)
        {
            var res = new Matrix(rows.Count, columns.Count);
            Parallel.For(0, rows.Count, i =>
            {
                for (int j = 0; j < columns.Count; j++)
                    res[i, j] = this[rows[i], columns[j]];
            });
            return res;
        }

        /// <summary>
        /// Returns Identity matrix with size n x n
        /// Added by MRB
        /// </summary>
        /// <param name="n"></param>
        /// <returns></returns>
        public static Matrix I(int n)
        {
            var Res = new Matrix(n, n);
            Parallel.For(0, n, i =>
            {
                //for (int j = 0; j < n; j++)
                Res[i, i] = 1; // i == j ? 1 : 0;
            });
            return Res;
        }

        /// <summary>
        /// Returns Zero matrix with size m x n
        /// Added by MRB
        /// </summary>
        /// <param name="n"></param>
        /// <returns></returns>
        public static Matrix Zeros(int m, int n)
        {
            var Res = new Matrix(m, n);
            /*Parallel.For(0, n, i =>
            {
                for (int j = 0; j < n; j++)
                    Res[i, j] = 0;
            });*/ // not needed
            return Res;
        }

        /// <summary>
        /// Added by MRB
        /// </summary>
        /// <param name="fileName"></param>
        /// <param name="varName"></param>
        /// <param name="append"></param>
        /// <returns></returns>
        public bool Export4Matlab(string fileName, string varName = null, bool append = false, string format = "")
        {
            try
            {
                if (fileName.Contains("/") || fileName.Contains("\\"))
                    try
                    {
                        Directory.CreateDirectory(Path.GetDirectoryName(fileName));
                    }
                    catch { }
                using (var A = new System.IO.StreamWriter(fileName, append))
                {
                    if (varName != null)
                        A.Write(varName + "=[");
                    for (int i = 0; i < this.RowCount; i++)
                    {
                        for (int j = 0; j < this.ColumnCount; j++)
                            A.Write(this[i, j].ToString(format) + "\t");
                        A.WriteLine(";");
                    }
                    if (varName != null)
                        A.WriteLine("];");
                    A.Close();
                }
                return true;
            }
            catch (Exception)
            {
                return false;
            }
        }
        /// <summary>
        /// added by mrb
        /// </summary>
        /// <param name="fileName"></param>
        /// <returns></returns>
        public static Matrix Import(string fileName)
        {
            var L = File.ReadAllLines(fileName);
            var A = new List<List<double>>();
            foreach (var l in L)
                if (l.Trim() != "")
                {
                    var V = l.Split(new char[] { '\t', ',', ' ' }, StringSplitOptions.RemoveEmptyEntries);
                    var a = new List<double>();
                    foreach (var v in V)
                        a.Add(Convert.ToDouble(v));
                    A.Add(a);
                }
            var res = new Matrix(A.Count, A[0].Count);
            for (int i = 0; i < res.RowCount; i++)
                for (int j = 0; j < res.ColumnCount; j++)
                    res[i, j] = A[i][j];
            return res;
        }

        public void Export(string fileName, string format = "")
        {
            if (fileName.Contains("/") || fileName.Contains("\\"))
                try
                {
                    Directory.CreateDirectory(Path.GetDirectoryName(fileName));
                }
                catch { }
            using (var f = new StreamWriter(fileName))
            {
                for (int i = 0; i < RowCount; i++)
                {
                    for (int j = 0; j < ColumnCount; j++)
                        f.Write(this[i, j].ToString() + (j < ColumnCount - 1 ? "\t" : ""));
                    if (i < RowCount - 1)
                        f.WriteLine();
                }
            }
        }
        /// <summary>
        /// Added by MRB
        /// </summary>
        /// <param name="fileName"></param>
        /// <param name="SaveLayout"></param>
        /// <returns></returns>
        public bool Export4Tecplot(string fileName, bool SaveLayout = true, bool OverrideLayout = false)
        {
            try
            {
                if (fileName.Contains("/") || fileName.Contains("\\"))
                    try
                    {
                        Directory.CreateDirectory(Path.GetDirectoryName(fileName));
                    }
                    catch { }
                using (var A = new System.IO.StreamWriter(fileName))
                {
                    A.WriteLine("VARIABLES = row col v");
                    for (int i = 0; i < this.RowCount; i++)
                    {
                        for (int j = 0; j < this.ColumnCount; j++)
                            A.WriteLine(i + "," + j + "," + this[i, j]);
                    }
                    A.Close();
                }
                if (SaveLayout && (OverrideLayout || !System.IO.File.Exists(System.IO.Path.ChangeExtension(fileName, "lay"))))
                {
                    var fileName_ = System.IO.Path.GetFileName(fileName);
                    #region Layout String
                    var layout_str = @"#!MC 1200
$!VarSet |LFDSFN1| = '" + fileName_ + @"'
$!VarSet |LFDSVL1| = 'row col v'
$!SETSTYLEBASE FACTORY
$!GLOBALLINKING 
  LINKCOLORMAPS = YES
$!PLOTOPTIONS 
  SUBDIVIDEALLCELLS = NO
$!GLOBALPAPER 
  PAPERSIZEINFO
    {
    LETTER
      {
      WIDTH = 8.5 HEIGHT = 11
      LEFTHARDCLIPOFFSET = 0.125
      RIGHTHARDCLIPOFFSET = 0.125
      TOPHARDCLIPOFFSET = 0.125
      BOTTOMHARDCLIPOFFSET = 0.125
      }
    }
$!PAGE 
  NAME = ''
  PAPERATTRIBUTES
    {
    BACKGROUNDCOLOR = WHITE
    ISTRANSPARENT = YES
    ORIENTPORTRAIT = NO
    SHOWGRID = YES
    SHOWRULER = NO
    SHOWPAPER = NO
    PAPERSIZE = LETTER
    RULERSPACING = ONEINCH
    PAPERGRIDSPACING = HALFINCH
    REGIONINWORKAREA { X1 = 1 Y1 = 0.25 X2 = 10 Y2 = 8.25 }
    }
### Frame Number 1 ###
$!READDATASET  '|LFDSFN1|'
  INITIALPLOTTYPE = CARTESIAN2D
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  ASSIGNSTRANDIDS = YES
  VARLOADMODE = BYNAME
  VARNAMELIST = '|LFDSVL1|'
$!REMOVEVAR |LFDSVL1|
$!REMOVEVAR |LFDSFN1|
$!FRAMELAYOUT 
  SHOWHEADER = NO
  HEADERCOLOR = RED
  XYPOS { X = 1 Y = 0.25 }
  WIDTH = 9
  HEIGHT = 8
$!THREEDAXIS 
  ASPECTRATIOLIMIT = 25
  BOXASPECTRATIOLIMIT = 25
$!PLOTTYPE  = CARTESIAN2D
$!FRAMENAME  = 'Frame 001'
$!ACTIVEFIELDMAPS  =  [1]
$!GLOBALRGB 
  RANGEMIN = 0
  RANGEMAX = 1
$!GLOBALCONTOUR  1
  VAR = 3
  LEGEND
    {
    SHOW = YES
    XYPOS { X = 98.5 }
    NUMBERTEXTSHAPE { HEIGHT = 1.5 }
    AUTORESIZE = YES
    LABELINCREMENT = 930.8
    }
  COLORCUTOFF { RANGEMIN = 0.25 RANGEMAX = 0.75 }
  COLORMAPFILTER { CONTINUOUSCOLOR { CMIN = 0 CMAX = 1 } }
$!FIELDMAP  [1]
  SCATTER { COLOR = MULTI FILLMODE = USELINECOLOR FRAMESIZE = " + (97.5 / (double)(Math.Max(RowCount, ColumnCount) - 1)) + @" }
  EDGELAYER { SHOW = NO }
  POINTS { POINTSTOPLOT = SURFACENODES }
  SURFACES { SURFACESTOPLOT = NONE }
  EFFECTS { LIGHTINGEFFECT = GOURAUD }
$!TWODAXIS 
  XDETAIL { VARNUM = 2 AXISLINE { AXISALIGNMENT = WITHGRIDMAX } }
  YDETAIL { VARNUM = 1 ISREVERSED = YES }
  VIEWPORTPOSITION { X1 = 9.5399 Y1 = 4.2886 Y2 = 91.893 }
  VIEWPORTTOPSNAPTARGET = 91.8926174497
$!VIEW FIT
$!FIELDLAYERS 
  SHOWMESH = NO
  SHOWSCATTER = YES
  SHOWSHADE = NO
  USETRANSLUCENCY = YES
$!STREAMTRACELAYERS 
  SHOW = NO
$!FRAMECONTROL ACTIVATEBYNUMBER
  FRAME = 1
$!SETSTYLEBASE CONFIG";
                    #endregion
                    using (var A = new System.IO.StreamWriter(System.IO.Path.ChangeExtension(fileName, "lay")))
                    {
                        A.WriteLine(layout_str);
                        A.Close();
                    }
                }
                return true;
            }
            catch (Exception)
            {
                return false;
            }
        }
        /// <summary>
        /// Added By MRB
        /// </summary>
        /// <param name="row_index"></param>
        /// <param name="row"></param>
        public void SetRow(int row_index, Vector row)
        {
            for (int i = 0; i < Math.Max(ColumnCount, row.Count); i++)
                this[row_index, i] = row[i];
        }
        /// <summary>
        /// Added By MRB
        /// </summary>
        /// <param name="col_index"></param>
        /// <param name="column"></param>
        public void SetColumn(int col_index, Vector column)
        {
            for (int i = 0; i < Math.Max(RowCount, column.Count); i++)
                this[i, col_index] = column[i];
        }


        public double Max()
        {
            if (_Data.Length == 0)
                return Double.NaN;
            var max = _Data[0];
            if (_Data.Length > 1000)
            {
                int N = 4;
                var max_ = new double[N];
                Parallel.For(0, N, i =>
                {
                    max_[i] = double.MinValue;
                    for (int j = i; j < _Data.Length; j += N)
                        max_[i] = Math.Max(max_[i], _Data[j]);
                });
                for (int i = 0; i < N; i++)
                    max = Math.Max(max, max_[i]);
            }
            else
                foreach (var d in Data)
                    max = Math.Max(max, d);
            return max;
        }
        public double MaxABS()
        {
            if (_Data.Length == 0)
                return Double.NaN;
            var max = 0.0;
            if (_Data.Length > 1000)
            {
                int N = 4;
                var max_ = new double[N];
                Parallel.For(0, N, i =>
                {
                    for (int j = i; j < _Data.Length; j += N)
                        max_[i] = Math.Max(Math.Abs(max_[i]), _Data[j]);
                });
                for (int i = 0; i < N; i++)
                    max = Math.Max(max, max_[i]);
            }
            else
                foreach (var d in Data)
                    max = Math.Max(Math.Abs(max), d);
            return max;
        }
        public double Min()
        {
            if (_Data.Length == 0)
                return Double.NaN;
            var min = _Data[0];
            if (_Data.Length > 1000)
            {
                int N = 4;
                var min_ = new double[N];
                Parallel.For(0, N, i =>
                {
                    min_[i] = double.MaxValue;
                    for (int j = i; j < _Data.Length; j += N)
                        min_[i] = Math.Min(min_[i], _Data[j]);
                });
                for (int i = 0; i < N; i++)
                    min = Math.Min(min, min_[i]);
            }
            else
                foreach (var d in Data)
                    min = Math.Min(min, d);
            return min;
        }
        public Matrix Abs()
        {
            var res = new Matrix(RowCount, ColumnCount);
            Parallel.For(0, _Data.Length, i =>
            {
                res._Data[i] = Math.Abs(_Data[i]);
            });
            return res;
        }

        public Matrix Round(int digits = 0)
        {
            var res = new Matrix(RowCount, ColumnCount);
            Parallel.For(0, _Data.Length, i =>
            {
                res._Data[i] = Math.Round(_Data[i], digits);
            });
            return res;
        }
        #endregion

        #region Matrix-Matrix Multiplication

        /// <summary>
        /// Matrix multiplication.
        /// </summary>
        /// <param name="A"> The left side matrix of the multiplication operator.</param>
        /// <param name="B">The right side matrix of the multiplication operator.</param>
        /// <returns>A matrix that represents the result of the matrix multiplication.</returns>
        public static Matrix operator *(BaseMatrix A, BaseMatrix B)
        {
            return A.Multiply(B);
        }

        #endregion

        #region Matrix-Matrix Addition

        /// <summary>
        /// Matrix addition.
        /// </summary>
        /// <param name="A">The left side matrix of the addition operator.</param>
        /// <param name="B">The right side matrix of the addition operator.</param>
        /// <returns>A matrix that represents the result of the matrix addition.</returns>
        public static Matrix operator +(BaseMatrix A, BaseMatrix B)
        {
            return A.Add(B);
        }


        #endregion

        #region Matrix-Matrix Subtraction

        /// <summary>
        /// Matrix subtraction.
        /// </summary>
        /// <param name="A"> The left side matrix of the subtraction operator.</param>
        /// <param name="B">The right side matrix of the subtraction operator.</param>
        /// <returns>A matrix that represents the result of the matrix subtraction.</returns>
        public static Matrix operator -(BaseMatrix A, BaseMatrix B)
        {
            return A.Subtract(B);
        }


        #endregion


        #region IMatrix<double> Members

        /// <summary>
        /// Copy all elements of this matrix to a rectangular 2D array.
        /// </summary>
        /// <returns>A rectangular 2D array.</returns>
        public double[,] CopyToArray()
        {
            double[,] matrixData = new double[this._RowCount, this._ColumnCount];

            for (int j = 0; j < this._ColumnCount; j++)
            {
                for (int i = 0; i < this._RowCount; i++)
                {
                    matrixData[i, j] = this._Data[i + j * this._RowCount];
                }
            }

            return matrixData;
        }

        /// <summary>
        /// Copy all elements of this matrix to a jagged array.
        /// </summary>
        /// <returns>A jagged array.</returns>
        public double[][] CopyToJaggedArray()
        {

            double[][] newData = new double[this._RowCount][];
            for (int i = 0; i < this._RowCount; i++)
            {
                double[] row = new double[this._ColumnCount];
                for (int j = 0; j < this._ColumnCount; j++)
                {
                    row[j] = this._Data[i + j * this._RowCount];
                }

                newData[i] = row;
            }

            return newData;
        }

        public ComplexMatrix CopyToComplex()
        {
            ComplexMatrix complexMatrix = new ComplexMatrix(this._RowCount, this._ColumnCount);
            Complex[] data = complexMatrix.Data;
            for (int i = 0; i < this._Data.Length; i++)
            {
                data[i].Real = this._Data[i];
            }

            return complexMatrix;
        }

        #endregion
    }
}
