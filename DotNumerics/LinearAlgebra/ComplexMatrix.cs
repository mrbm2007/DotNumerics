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
using DotNumerics.FortranLibrary;
using System.Threading.Tasks;
using System.Linq;

namespace DotNumerics.LinearAlgebra
{
    //[DebuggerDisplay("[{RowCount},{ColumnCount}]", Name = "MatrixComplex")]
    /// <summary>
    /// Represents a Complex Matrix.
    /// </summary>
    [DebuggerDisplay("[{RowCount},{ColumnCount}]")]
    [DebuggerTypeProxy(typeof(MatrixComplexDebuggerDisplay))]
    public class ComplexMatrix : IMatrix<Complex>
    {
        #region Fields
        /// <summary>
        /// Los datos de la matriz, los datos se almacenan en un un array unidimensional,
        /// Los elementos se almacenan por columnas, esto para que sean compatible con los Arrays de Fortran
        /// </summary>
        [DebuggerBrowsable(DebuggerBrowsableState.Never)]
        protected Complex[] _Data;
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
        /// Initializes a new instance of the MatrixComplex class of the given size.
        /// </summary>
        /// <param name="rows">Number of rows.</param>
        /// <param name="columns">Number of columns.</param>
        public ComplexMatrix(int rows, int columns)
        {
            if (!BaseMatrix.AllowEmptyMatrices)
            {
                if (rows < 1) throw new System.ArgumentException("rows < 1");
                if (columns < 1) throw new System.ArgumentException("columns < 1");
            }

            this._Data = new Complex[rows * columns];
            this._RowCount = rows;
            this._ColumnCount = columns;
            for (int i = 0; i < this._Data.Length; i++)
                this._Data[i] = new Complex(0, 0);
        }

        /// <summary>
        /// Initializes a new instance of the MatrixComplex class of the given size using a array
        /// </summary>
        /// <param name="rows">Number of rows.</param>
        /// <param name="columns">Number of columns.</param>
        /// <param name="Data">The data</param>
        internal ComplexMatrix(int rows, int columns, Complex[] Data)
        {
            if (!BaseMatrix.AllowEmptyMatrices)
            {
                if (rows < 1) throw new System.ArgumentException("rows < 1");
                if (columns < 1) throw new System.ArgumentException("columns < 1");
            }
            this._Data = new Complex[rows * columns];
            this._RowCount = rows;
            this._ColumnCount = columns;
            //Si incluye la posibilidad de que los datos tengan menos valores que la matriz a crear 
            for (int i = 0; i < Math.Min(this._Data.Length, Data.Length); i++)
            {
                this._Data[i] = Data[i];
            }
        }

        /// <summary>
        /// Initializes a new instance of the MatrixComplex class of the given size.
        /// </summary>
        /// <param name="size">Size</param>
        public ComplexMatrix(int size)
        {
            if (!BaseMatrix.AllowEmptyMatrices)
                if (size < 1) throw new System.ArgumentException("size < 1");

            this._Data = new Complex[size * size];
            this._RowCount = size;
            this._ColumnCount = size;
            for (int i = 0; i < this._Data.Length; i++)
                this._Data[i] = new Complex(0, 0);
        }

        /// <summary>
        /// Initializes a new instance of the MatrixComplex class of the given size using a array
        /// </summary>
        /// <param name="size">Size</param>
        /// <param name="Data">The data</param>
        internal ComplexMatrix(int size, Complex[] Data)
        {
            if (!BaseMatrix.AllowEmptyMatrices)
                if (size < 1) throw new System.ArgumentException("size < 1");

            this._Data = new Complex[size * size];
            this._RowCount = size;
            this._ColumnCount = size;
            //Si incluye la posibilidad de que los datos tengan menos valores que la matriz a crear 
            for (int i = 0; i < Math.Min(this._Data.Length, Data.Length); i++)
            {
                this._Data[i] = Data[i];
            }
        }

        public ComplexMatrix(Matrix real, Matrix imag) :
            this(real.RowCount, real.ColumnCount, real.Data, imag.Data)
        {
            if (real.RowCount != imag.RowCount)
                throw new ArgumentException("ComplexMatrix constructor, real and imag RowCount not equal: " + real.RowCount + "," + imag.RowCount);
        }
        public ComplexMatrix(int rows, int columns, double[] real, double[] imag)
        {
            if (!BaseMatrix.AllowEmptyMatrices)
            {
                if (rows < 1) throw new ArgumentException("rows < 1");
                if (columns < 1) throw new ArgumentException("columns < 1");
            }
            if (real.Length != imag.Length)
                throw new ArgumentException("ComplexMatrix constructor, real and imag length not equal: " + real.Length + "," + imag.Length);
            this._RowCount = rows;
            this._ColumnCount = columns;
            this._Data = new Complex[rows * columns];
            for (int i = 0; i < this._Data.Length; i++)
            {
                this._Data[i].Real = real[i];
                this._Data[i].Imaginary = imag[i];
            }
        }
        public ComplexMatrix(int rows, int columns, Func<int, int, Complex> filler) : this(rows, columns)
        {
            Fill(filler);
        }
        public ComplexMatrix(int rows, int columns,
            Func<int, int, double> filler_r, Func<int, int, double> filler_i) :
            this(rows, columns)
        {
            Fill(filler_r, filler_i);
        }
        public ComplexMatrix(int rows, int columns, Func<int, Complex> filler) : this(rows, columns)
        {
            Fill_Diag(filler);
        }

        public void Fill(Func<int, int, Complex> filler)
        {
            Parallel.For(0, ColumnCount, j =>
            {
                var jj = j * this._RowCount;
                for (int i = 0; i < RowCount; i++)
                    Data[jj + i] = filler(i, j);
            });
        }
        public void Fill(Func<int, int, double> filler_r, Func<int, int, double> filler_i)
        {
            Parallel.For(0, ColumnCount, j =>
            {
                var jj = j * this._RowCount;
                for (int i = 0; i < RowCount; i++)
                {
                    Data[jj + i].Real = filler_r(i, j);
                    Data[jj + i].Imaginary = filler_i(i, j);
                }
            });
        }
        public void Fill_Diag(Func<int, Complex> filler)
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
        public Complex[] Data
        {
            get { return this._Data; }
        }

        public double[] Data_r
        {
            get
            {
                return Data.Select(d => d.Real).ToArray();
            }
            set
            {
                for (int i = 0; i < Data.Length; i++)
                    Data[i].Real = value[i];
            }
        }
        public double[] Data_i
        {
            get
            {
                return Data.Select(d => d.Imaginary).ToArray();
            }
            set
            {
                for (int i = 0; i < Data.Length; i++)
                    Data[i].Imaginary = value[i];
            }
        }

        /// <summary>
        /// Returns the number of rows.
        /// </summary>
        [DebuggerBrowsable(DebuggerBrowsableState.Never)]
        public int RowCount
        {
            get { return _RowCount; }
            set { _RowCount = value; }
        }
        /// <summary>
        /// Returns the number of columns.
        /// </summary>
        [DebuggerBrowsable(DebuggerBrowsableState.Never)]
        public int ColumnCount
        {
            get { return _ColumnCount; }
            set { _ColumnCount = value; }
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
        /// Returns the value of a element of the matrix.
        /// </summary>
        /// <param name="row">The row value (zero-based).</param>
        /// <param name="column">The column value (zero-based).</param>
        /// <returns>The matrix value at (row, column).</returns>
        public virtual Complex this[int row, int column]
        {
            get
            {
                if (column >= this._ColumnCount)
                {
                    throw new ArgumentException("Index was outside the bounds of the matrix.");
                }
                return this._Data[row + column * this._RowCount];
            }
            set
            {
                if (column >= this._ColumnCount)
                {
                    throw new ArgumentException("Index was outside the bounds of the matrix.");
                }
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
        public ComplexMatrix this[int startRow, int startCol, int rowCount, int colCount]
        {
            get { return SubMatrix(startRow, startCol, rowCount, colCount); }
            set { InsertAt(startRow, startCol, value, rowCount, colCount); }
        }
        public ComplexMatrix this[List<int> Rows, List<int> Columns]
        {
            get { return SubMatrix(Rows, Columns); }
            set { InsertAt(Rows, Columns, value); }
        }

        public RealProp real => new RealProp(this);
        public ImagProp imag => new ImagProp(this);

        public class RealProp
        {
            ComplexMatrix owner;
            public RealProp(ComplexMatrix owner)
            {
                this.owner = owner;
            }
            public virtual double this[int row, int column]
            {
                get { return owner._Data[row + column * owner._RowCount].Real; }
                set { owner._Data[row + column * owner._RowCount].Real = value; }
            }
        }
        public class ImagProp
        {
            ComplexMatrix owner;
            public ImagProp(ComplexMatrix owner)
            {
                this.owner = owner;
            }
            public virtual double this[int row, int column]
            {
                get { return owner._Data[row + column * owner._RowCount].Imaginary; }
                set { owner._Data[row + column * owner._RowCount].Imaginary = value; }
            }
        }

        #endregion

        #region	 Private Methods

        /// <summary>Check if size(this) == size(B) </summary>
        private void CheckMatrixDimensions(ComplexMatrix B)
        {
            if (this._RowCount != B.RowCount || B.ColumnCount != this._ColumnCount)
            {
                throw new System.ArgumentException("Matrix dimensions must agree.");
            }
        }

        /// <summary>Check if size(this) == size(B) </summary>
        private void CheckMatrixDimensions(Matrix B)
        {
            if (this._RowCount != B.RowCount || B.ColumnCount != this._ColumnCount)
            {
                throw new System.ArgumentException("Matrix dimensions must agree.");
            }
        }
        #endregion //  Private Methods

        #region Elementary linear operations

        ///// <summary>
        ///// aij=Math.Abs(aij)
        ///// </summary>
        //public virtual void ElementsAbs()
        //{
        //    for (int i = 0; i < this.MeData.Length; i++)
        //    {
        //        this.MeData[i] =Complex.  Math.Abs(this.MeData[i]);
        //    }
        //}

        ///// <summary>
        ///// Element-by-element division: aij = aij/bij
        ///// </summary>
        ///// <param name="B"></param>
        //public virtual void ElemntsDiv(MatrixComplex B)
        //{
        //    CheckMatrixDimensions(B);
        //    Complex[] BData = B.Data;
        //    for (int i = 0; i < this.MeData.Length; i++)
        //    {
        //        this.MeData[i] /= BData[i];
        //    }
        //}

        /// <summary>
        /// Element-by-element multiplication: aij = aij*bij
        /// </summary>
        /// <param name="B">The B MatrixComplex</param>
        public virtual void ElemntsMult(ComplexMatrix B)
        {
            CheckMatrixDimensions(B);
            Complex[] BData = B.Data;
            for (int i = 0; i < this._Data.Length; i++)
            {
                this._Data[i] = this._Data[i] * BData[i];
            }
        }

        /// <summary>
        /// In place addition A=A+B
        /// </summary>
        /// <param name="B">The B MatrixComplex</param>
        public virtual void Add(ComplexMatrix B)
        {
            CheckMatrixDimensions(B);
            Complex[] BData = B.Data;
            for (int i = 0; i < this._Data.Length; i++)
            {
                this._Data[i] = this._Data[i] + BData[i];
            }
        }

        /// <summary>
        /// In place scalar-matrix multiplication, A=s*A
        /// </summary>
        /// <param name="s">The scalar s.</param>
        public virtual void Multiply(double s)
        {
            for (int i = 0; i < this._Data.Length; i++)
            {
                this._Data[i].Real = s * this._Data[i].Real;
                this._Data[i].Imaginary = s * this._Data[i].Imaginary;
            }
        }
        /// <summary>
        /// In place scalar-matrix multiplication, A=c*A
        /// </summary>
        /// <param name="s">The scalar s.</param>
        public virtual void MultiplyC(Complex c)
        {
            for (int i = 0; i < this._Data.Length; i++)
            {
                this._Data[i] = c * this._Data[i];
            }
        }


        /// <summary>
        /// In place matrix subtraction, A=A-B.
        /// </summary>
        /// <param name="B">The B MatrixComplex.</param>
        public virtual void Subtract(ComplexMatrix B)
        {
            CheckMatrixDimensions(B);
            Complex[] BData = B.Data;
            for (int i = 0; i < this._Data.Length; i++)
            {
                this._Data[i].Real = this._Data[i].Real - BData[i].Real;
                this._Data[i].Imaginary = this._Data[i].Imaginary - BData[i].Imaginary;
            }
        }

        /// <summary>
        /// In place unary minus -A.
        /// </summary>
        public virtual void UnaryMinus()
        {
            for (int i = 0; i < this._Data.Length; i++)
            {
                this._Data[i].Real = -this._Data[i].Real;
                this._Data[i].Imaginary = -this._Data[i].Imaginary;
            }
        }

        #endregion

        #region Methods


        /// <summary>
        /// Added By MRB
        /// </summary>
        /// <param name="startRow"></param>
        /// <param name="startCol"></param>
        /// <param name="M"></param>
        /// <param name="rowCount"></param>
        /// <param name="colCount"></param>
        public void InsertAt(int startRow, int startCol, ComplexMatrix M, int rowCount = -1, int colCount = -1)
        {
            if (rowCount == -1)
                rowCount = M.RowCount;
            if (colCount == -1)
                colCount = M.ColumnCount;
            for (int i = 0; i < rowCount; i++)
                for (int j = 0; j < colCount; j++)
                    this[i + startRow, j + startCol] = M[i, j];
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
            for (int i = 0; i < rowCount; i++)
                for (int j = 0; j < colCount; j++)
                    this[i + startRow, j + startCol] = M[i, j];
        }

        public void InsertAt(List<int> Rows, List<int> Columns, ComplexMatrix A1)
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
        public void PlusAt(int startRow, int startCol, ComplexMatrix M, int rowCount = -1, int colCount = -1)
        {
            if (rowCount == -1)
                rowCount = M.RowCount;
            if (colCount == -1)
                colCount = M.ColumnCount;
            for (int i = 0; i < rowCount; i++)
                for (int j = 0; j < colCount; j++)
                    this[i + startRow, j + startCol] += M[i, j];
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
            for (int i = 0; i < rowCount; i++)
                for (int j = 0; j < colCount; j++)
                    this[i + startRow, j + startCol] += M[i, j];
        }



        /// <summary>
        ///  Gets a column vector of this matrix at the selected position.
        /// </summary>
        /// <param name="columnIndex">The column index (zero-based).</param>
        /// <returns>The column vector.</returns>
        public ComplexVector GetColumnVector(int columnIndex)
        {

            if (columnIndex >= this.ColumnCount)
                throw new System.ArgumentException("columnIndex >= number of columns.");

            if (columnIndex < 0)
                throw new System.ArgumentException("columnIndex < 0");


            var columnVect = new ComplexVector(VectorType.Column, this._RowCount);
            Array.Copy(this._Data, columnIndex * this._RowCount, columnVect.Data, 0, this._RowCount);
            return columnVect;
        }

        public ComplexMatrix SubMatrix(int startRow = 0, int startCol = 0, int rowCount = -1, int colCount = -1)
        {
            if (rowCount == -1)
                rowCount = RowCount - startRow;
            if (colCount == -1)
                colCount = ColumnCount - startCol;
            var res = new ComplexMatrix(rowCount, colCount);
            for (int i = 0; i < rowCount; i++)
                for (int j = 0; j < colCount; j++)
                    res[i, j] = this[i + startRow, j + startCol];
            return res;
        }

        public ComplexMatrix SubMatrix(List<int> rows, List<int> columns)
        {
            var res = new ComplexMatrix(rows.Count, columns.Count);
            Parallel.For(0, rows.Count, i =>
            {
                for (int j = 0; j < columns.Count; j++)
                    res[i, j] = this[rows[i], columns[j]];
            });
            return res;
        }

        public void Resize(int newRowCount, int newColumnCount, bool keep_elems = true)
        {
            if (keep_elems && newRowCount != this.RowCount)
            {
                var tmp = new ComplexMatrix(this._RowCount, this._ColumnCount);
                Array.Copy(_Data, tmp._Data, _Data.Length);
                _Data = new Complex[newRowCount * newColumnCount];
                this._RowCount = newRowCount;
                this._ColumnCount = newColumnCount;

                System.Threading.Tasks.Parallel.For(0, Math.Min(tmp.RowCount, this.RowCount), i =>
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
                var A = new System.IO.StreamWriter(fileName, append);
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
                return true;
            }
            catch (Exception)
            {
                return false;
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
                using (var A = new System.IO.StreamWriter(fileName))
                {
                    A.WriteLine("VARIABLES = row col real imag");
                    for (int i = 0; i < this.RowCount; i++)
                    {
                        for (int j = 0; j < this.ColumnCount; j++)
                            A.WriteLine(i + "," + j + "," + this[i, j].Real + "," + this[i, j].Imaginary);
                    }
                    A.Close();
                }
                if (SaveLayout && (OverrideLayout || !System.IO.File.Exists(System.IO.Path.ChangeExtension(fileName, "lay"))))
                {
                    var fileName_ = System.IO.Path.GetFileName(fileName);
                    #region Layout String
                    var layout_str = @"#!MC 1200
$!VarSet |LFDSFN1| = '" + fileName_ + @"'
$!VarSet |LFDSVL1| = 'row col real imag'
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
        /// Gets the column vectors of this matrix.
        /// </summary>
        /// <returns>The columns vectors.</returns>
        public ComplexVector[] GetColumnVectors()
        {
            var res = new ComplexVector[ColumnCount];
            for (int i = 0; i < ColumnCount; i++)
            {
                var columnVect = new ComplexVector(VectorType.Column, this._RowCount);
                Array.Copy(this._Data, i * this._RowCount, columnVect.Data, 0, this._RowCount);
                res[i] = columnVect;
            }
            return res;
        }

        /// <summary>
        /// Gets the row vectors of this matrix.
        /// </summary>
        /// <returns>The row vectors.</returns>
        public ComplexVector[] GetRowVectors()
        {
            ComplexVector[] rowVects = new ComplexVector[this.RowCount];

            Complex[] VectData;
            for (int i = 0; i < this._RowCount; i++)
            {
                rowVects[i] = new ComplexVector(VectorType.Row, this._ColumnCount);
                VectData = rowVects[i].Data;
                for (int j = 0; j < VectData.Length; j++)
                {
                    VectData[j] = this._Data[i + j * this._RowCount];
                }
            }

            return rowVects;
        }

        /// <summary>
        /// Gets a matrix that contains the real part of this matrix.
        /// </summary>
        /// <returns>A matrix that contains the real part of this matrix. </returns>
        public Matrix GetReal()
        {
            Matrix RealMatrix = new Matrix(this.RowCount, this.ColumnCount);
            double[] RealData = RealMatrix.Data;
            for (int i = 0; i < this._Data.Length; i++)
            {
                RealData[i] = this._Data[i].Real;
            }
            return RealMatrix;
        }

        /// <summary>
        /// Gets a matrix that contains the imaginary part of this matrix.
        /// </summary>
        /// <returns>A matrix that contains the imaginary part of this matrix. </returns>
        public Matrix GetImag()
        {
            Matrix ImagMatrix = new Matrix(this.RowCount, this.ColumnCount);
            double[] ImagData = ImagMatrix.Data;
            for (int i = 0; i < this._Data.Length; i++)
            {
                ImagData[i] = this._Data[i].Imaginary;
            }
            return ImagMatrix;
        }


        /// <summary>
        /// Sets the real part of the elements of this matrix equal to the elemnets of a real matrix.
        /// </summary>
        /// <param name="RM">A matrix that contains the values of the real part.</param>
        public void SetReal(Matrix RM)
        {
            this.CheckMatrixDimensions(RM);
            double[] RealData = RM.Data;
            for (int i = 0; i < this._Data.Length; i++)
            {
                this._Data[i].Real = RealData[i];
            }
        }

        /// <summary>
        /// Sets the imaginary part of the elements of this matrix equal to the elemnets of a real matrix.
        /// </summary>
        /// <param name="IM">A matrix that contains the values of the imaginary part.</param>
        public void SetImag(Matrix IM)
        {
            this.CheckMatrixDimensions(IM);
            double[] ImagData = IM.Data;
            for (int i = 0; i < this._Data.Length; i++)
            {
                this._Data[i].Imaginary = ImagData[i];
            }
        }


        /// <summary>
        /// Returns the string of the  matrix.
        /// </summary>
        /// <returns>The string of the  matrix.</returns>
        public string MatrixToString()
        {
            using (StringWriter writer = new StringWriter())
            {
                for (int i = 0; i < this._RowCount; i++)
                {
                    for (int j = 0; j < this._ColumnCount; j++)
                        writer.Write(this[i, j] + ", ");
                    writer.WriteLine();
                }
                return writer.ToString();
            }
        }
        /// <summary>
        /// Returns the string of the  matrix.
        /// </summary>
        /// <param name="format">A numeric format string.</param>
        /// <returns>The string of the  matrix.</returns>
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


        ///// <summary>
        ///// maximum column sum.
        ///// </summary>
        ///// <returns>maximum column sum.</returns>
        //public double Norm1()
        //{
        //    double n = 0.0;
        //    double ColSum = 0.0;
        //    int NRows = this.MeRowCount;

        //    for (int j = 0; j < this.MeColumnCount; j++)
        //    {
        //        ColSum = 0.0;
        //        for (int i = 0; i < this.MeRowCount; i++)
        //        {
        //            ColSum += Math.Abs(this.MeData[i + j * NRows]);

        //        }
        //        n = Math.Max(n, ColSum);
        //    }
        //    return n;
        //}
        ///// <summary>
        ///// 
        ///// </summary>
        ///// <returns></returns>
        //public double InfinityNorm()
        //{
        //    double n = 0.0;
        //    double RowSum = 0.0;
        //    int NRows = this.MeRowCount;
        //    for (int i = 0; i < this.MeRowCount; i++)
        //    {
        //        RowSum = 0.0;
        //        for (int j = 0; j < this.MeColumnCount; j++)
        //        {
        //            RowSum += Math.Abs(this.MeData[i + j * NRows]);
        //        }
        //        n = Math.Max(n, RowSum);
        //    }
        //    return n;
        //}

        ///// <summary>Frobenius norm</summary>
        ///// <returns>Sqrt of sum of squares of all elements.</returns>
        //public double FrobeniusNorm()
        //{
        //    double n=0;
        //    for(int i=0; i<this.MeData.Length; i++)
        //    {
        //        n=this.Hypot(n,this.MeData[i]);
        //    }
        //    return n;
        //}

        ///// <summary>sqrt(a^2 + b^2) without under/overflow.</summary>
        //private  double Hypot(double a, double b) 
        //{
        //    double r;
        //    if (Math.Abs(a) > Math.Abs(b)) 
        //    {
        //        r = b/a;
        //        r = Math.Abs(a) * Math.Sqrt(1 + r * r);
        //    } 
        //    else if (b != 0) 
        //    {
        //        r = a/b;
        //        r = Math.Abs(b) * Math.Sqrt(1 + r * r);
        //    } 
        //    else 
        //    {
        //        r = 0.0;
        //    }
        //    return r;
        //}


        #endregion

        #region Matrix-Matrix Multiplication

        /// <summary>
        /// Matrix multiplication.
        /// </summary>
        /// <param name="A"> The left side matrix of the multiplication operator.</param>
        /// <param name="B">The right side matrix of the multiplication operator.</param>
        /// <returns>A matrix that represents the result of the matrix multiplication.</returns>
        public static ComplexMatrix operator *(ComplexMatrix A, ComplexMatrix B)
        {
            if (B.RowCount != A.ColumnCount)
            {
                throw new System.ArgumentException("Matrix dimensions are not valid.");
            }

            if (ARMADILO.Enabled)
                return ARMADILO.mat_multiply(A, B);

            ComplexMatrix C = new ComplexMatrix(A.RowCount, B.ColumnCount);

            Complex[] AData = A.Data;
            Complex[] BData = B.Data;
            Complex[] CData = C.Data;

            int ARows = A.RowCount;
            int AColumns = A.ColumnCount;

            int BRows = B.RowCount;
            int BColumns = B.ColumnCount;

            for (int j = 0; j < BColumns; j++)
            {
                for (int i = 0; i < ARows; i++)
                {
                    //Sum.Imaginary = 0.0;
                    //Sum.Real = 0.0;
                    Complex Sum = new Complex(0.0, 0.0);
                    for (int k = 0; k < AColumns; k++)
                    {
                        Sum += AData[i + k * ARows] * BData[k + j * BRows];
                    }
                    CData[i + j * ARows] = Sum;
                }
            }
            return C;
        }
        /// <summary>
        /// Matrix multiplication.
        /// </summary> 
        public static ComplexMatrix operator *(ComplexMatrix A, BaseMatrix B)
        {
            if (B.RowCount != A.ColumnCount)
            {
                throw new System.ArgumentException("Matrix dimensions are not valid.");
            }
            if (ARMADILO.Enabled && B is Matrix)
                return ARMADILO.mat_multiply(A, B as Matrix);

            ComplexMatrix C = new ComplexMatrix(A.RowCount, B.ColumnCount);

            Complex[] AData = A.Data;
            double[] BData = B.Data;
            Complex[] CData = C.Data;

            int ARows = A.RowCount;
            int AColumns = A.ColumnCount;

            int BRows = B.RowCount;
            int BColumns = B.ColumnCount;

            for (int j = 0; j < BColumns; j++)
            {
                for (int i = 0; i < ARows; i++)
                {
                    //Sum.Imaginary = 0.0;
                    //Sum.Real = 0.0;
                    Complex Sum = new Complex(0.0, 0.0);
                    for (int k = 0; k < AColumns; k++)
                    {
                        Sum += AData[i + k * ARows] * BData[k + j * BRows];
                    }
                    CData[i + j * ARows] = Sum;
                }
            }
            return C;
        }

        /// <summary>
        /// Matrix multiplication.
        /// </summary>
        /// <param name="A"> The left side matrix of the multiplication operator.</param>
        /// <param name="B">The right side matrix of the multiplication operator.</param>
        /// <returns>A matrix that represents the result of the matrix multiplication.</returns>
        public static ComplexMatrix operator *(BaseMatrix A, ComplexMatrix B)
        {
            if (B.RowCount != A.ColumnCount)
            {
                throw new System.ArgumentException("Matrix dimensions are not valid.");
            }
            if (ARMADILO.Enabled && A is Matrix)
                return ARMADILO.mat_multiply(A as Matrix, B);

            ComplexMatrix C = new ComplexMatrix(A.RowCount, B.ColumnCount);

            double[] AData = A.Data;
            Complex[] BData = B.Data;
            Complex[] CData = C.Data;

            int ARows = A.RowCount;
            int AColumns = A.ColumnCount;

            int BRows = B.RowCount;
            int BColumns = B.ColumnCount;

            for (int j = 0; j < BColumns; j++)
            {
                for (int i = 0; i < ARows; i++)
                {
                    //Sum.Imaginary = 0.0;
                    //Sum.Real = 0.0;
                    Complex Sum = new Complex(0.0, 0.0);
                    for (int k = 0; k < AColumns; k++)
                    {
                        Sum += AData[i + k * ARows] * BData[k + j * BRows];
                    }
                    CData[i + j * ARows] = Sum;
                }
            }
            return C;
        }

        /// <summary>complex-Matrix multiplication.</summary>
        /// <param name="c"> The left side scalar of the multiplication operator.</param>
        /// <param name="B">The right side matrix of the multiplication operator.</param>
        /// <returns>A matrix that represents the result of the multiplication.</returns>
        public static ComplexMatrix operator *(Complex c, ComplexMatrix B)
        {

            ComplexMatrix C = new ComplexMatrix(B.RowCount, B.ColumnCount);

            Complex[] BData = B.Data;
            Complex[] CData = C.Data;


            for (int i = 0; i < BData.Length; i++)
            {
                CData[i] = c * BData[i];
            }
            return C;
        }

        public static ComplexVector operator *(ComplexMatrix A, Vector b)
        {
            if (b.Count != A.ColumnCount)
            {
                throw new System.ArgumentException("Matrix dimensions are not valid:" + A.RowCount + "x" + A.ColumnCount + " , " + b.Count);
            }
            if (ARMADILO.Enabled)
                return ARMADILO.mv_multiply(A, b);

            var C = new ComplexVector(A.RowCount);

            var AData = A.Data;
            var BData = b.Data;
            var CData = C.Data;

            int ARows = A.RowCount;
            int AColumns = A.ColumnCount;

            for (int i = 0; i < ARows; i++)
            {
                //Sum.Imaginary = 0.0;
                //Sum.Real = 0.0;
                Complex Sum = new Complex(0.0, 0.0);
                for (int k = 0; k < AColumns; k++)
                {
                    Sum += AData[i + k * ARows] * BData[k];
                }
                CData[i] = Sum;
            }
            return C;
        }
        #endregion

        #region Matrix-Matrix Addition
        /// <summary>
        /// Matrix addition.
        /// </summary>
        /// <param name="A">The left side matrix of the addition operator.</param>
        /// <param name="B">The right side matrix of the addition operator.</param>
        /// <returns>A matrix that represents the result of the matrix addition.</returns>
        public static ComplexMatrix operator +(ComplexMatrix A, ComplexMatrix B)
        {
            if (B.RowCount != A.RowCount || B.ColumnCount != A.ColumnCount)
            {
                throw new System.ArgumentException("Matrix dimensions are not valid.");
            }

            ComplexMatrix C = new ComplexMatrix(A.RowCount, A.ColumnCount);

            Complex[] AData = A.Data;
            Complex[] BData = B.Data;
            Complex[] CData = C.Data;

            for (int i = 0; i < AData.Length; i++)
            {
                CData[i] = AData[i] + BData[i];
            }

            return C;
        }
        public static ComplexMatrix operator +(ComplexMatrix A, BaseMatrix B)
        {
            if (B.RowCount != A.RowCount || B.ColumnCount != A.ColumnCount)
            {
                throw new System.ArgumentException("Matrix dimensions are not valid.");
            }

            ComplexMatrix C = new ComplexMatrix(A.RowCount, A.ColumnCount);

            var AData = A.Data;
            var BData = B.Data;
            var CData = C.Data;

            for (int i = 0; i < AData.Length; i++)
                CData[i] = AData[i] + BData[i];

            return C;
        }
        public static ComplexMatrix operator +(BaseMatrix A, ComplexMatrix B)
        {
            if (B.RowCount != A.RowCount || B.ColumnCount != A.ColumnCount)
            {
                throw new System.ArgumentException("Matrix dimensions are not valid.");
            }

            ComplexMatrix C = new ComplexMatrix(A.RowCount, A.ColumnCount);

            var AData = A.Data;
            var BData = B.Data;
            var CData = C.Data;

            for (int i = 0; i < AData.Length; i++)
                CData[i] = AData[i] + BData[i];

            return C;
        }


        /// <summary>
        /// Unary minus.
        /// </summary>
        /// <param name="A"> The Matrix</param>
        /// <returns>Matrix r[i] = -this[i]</returns>
        public static ComplexMatrix operator -(ComplexMatrix A)
        {
            return A.UnaryMinus_();
        }
        public ComplexMatrix UnaryMinus_()
        {
            var C = new ComplexMatrix(this._RowCount, this._ColumnCount);
            var dataC = C.Data;

            System.Threading.Tasks.Parallel.For(0, this._Data.Length, i =>
            {
                dataC[i] = -this._Data[i];
            });
            return C;
        }
        #endregion

        #region Matrix-Matrix Subtraction

        ///// <summary>Matrix Subtraction</summary>
        /// <summary>
        /// Matrix subtraction.
        /// </summary>
        /// <param name="A"> The left side matrix of the subtraction operator.</param>
        /// <param name="B">The right side matrix of the subtraction operator.</param>
        /// <returns>A matrix that represents the result of the matrix subtraction.</returns>
        public static ComplexMatrix operator -(ComplexMatrix A, ComplexMatrix B)
        {
            if (B.RowCount != A.RowCount || B.ColumnCount != A.ColumnCount)
            {
                throw new System.ArgumentException("Matrix dimensions are not valid.");
            }

            ComplexMatrix C = new ComplexMatrix(A.RowCount, A.ColumnCount);

            Complex[] AData = A.Data;
            Complex[] BData = B.Data;
            Complex[] CData = C.Data;

            for (int i = 0; i < AData.Length; i++)
            {
                CData[i] = AData[i] - BData[i];
            }

            return C;
        }
        public static ComplexMatrix operator -(ComplexMatrix A, BaseMatrix B)
        {
            if (B.RowCount != A.RowCount || B.ColumnCount != A.ColumnCount)
            {
                throw new System.ArgumentException("Matrix dimensions are not valid.");
            }

            ComplexMatrix C = new ComplexMatrix(A.RowCount, A.ColumnCount);

            var AData = A.Data;
            var BData = B.Data;
            var CData = C.Data;

            for (int i = 0; i < AData.Length; i++)
            {
                CData[i] = AData[i] - BData[i];
            }

            return C;
        }
        public static ComplexMatrix operator -(BaseMatrix A, ComplexMatrix B)
        {
            if (B.RowCount != A.RowCount || B.ColumnCount != A.ColumnCount)
            {
                throw new System.ArgumentException("Matrix dimensions are not valid.");
            }

            ComplexMatrix C = new ComplexMatrix(A.RowCount, A.ColumnCount);

            var AData = A.Data;
            var BData = B.Data;
            var CData = C.Data;

            for (int i = 0; i < AData.Length; i++)
            {
                CData[i] = AData[i] - BData[i];
            }

            return C;
        }

        #endregion

        /// <summary>
        /// Transposed matrix.
        /// </summary>
        /// <returns>The transposed matrix.</returns>
        public ComplexMatrix Transpose()
        {
            ComplexMatrix AT = new ComplexMatrix(this._ColumnCount, this._RowCount);
            int ATRows = AT.RowCount;
            int ATColumns = AT.ColumnCount;
            var ATData = AT.Data;
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

        #region IMatrix<Complex> Members

        /// <summary>
        /// Copy all elements of this matrix to a rectangular 2D array.
        /// </summary>
        /// <returns>A rectangular 2D array.</</returns>
        public Complex[,] CopyToArray()
        {
            Complex[,] matrixData = new Complex[this._RowCount, this._ColumnCount];

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
        public Complex[][] CopyToJaggedArray()
        {
            Complex[][] newData = new Complex[this._RowCount][];
            for (int i = 0; i < this._RowCount; i++)
            {
                Complex[] row = new Complex[this._ColumnCount];
                for (int j = 0; j < this._ColumnCount; j++)
                {
                    row[j] = this._Data[i + j * this._RowCount];
                }

                newData[i] = row;
            }

            return newData;
        }

        #endregion

        internal void SwapColumns(int i, int j)
        {
            var i1 = _RowCount * i;
            var j1 = _RowCount * j;
            for (int c = 0; c < _RowCount; c++)
            {
                var tmp = _Data[i1 + c];
                _Data[i1 + c] = _Data[j1 + c];
                _Data[j1 + c] = tmp;
            }
            /*Export4Matlab("e:\\k1.txt", "k");
            for (int r = 0; r < RowCount; r++)
            {
                var tmp = this[r, i];
                this[r, i] = this[r, j];
                this[r, j] = tmp;
            }
            Export4Matlab("e:\\k2.txt", "k");*/
        }

        /*internal ComplexMatrix Transpose()
        {
            var AT = new ComplexMatrix(this._ColumnCount, this._RowCount);
            int ATRows = AT.RowCount;
            int ATColumns = AT.ColumnCount;
            var ATData = AT._Data;
            for (int j = 0; j < this._ColumnCount; j++)
            {
                for (int i = 0; i < this._RowCount; i++)
                {
                    ATData[j + i * ATRows] = this._Data[i + j * this._RowCount];
                }
            }
            return AT;
        }*/

        public static ComplexMatrix operator /(ComplexMatrix A, double s)
        {
            var C = new ComplexMatrix(A._RowCount, A._ColumnCount);

            var dataC = C.Data;
            for (int i = 0; i < A._Data.Length; i++)
            {
                dataC[i] = A._Data[i] / s;
            }
            return C;
        }
        public static ComplexMatrix operator /(ComplexMatrix A, Complex s)
        {
            var C = new ComplexMatrix(A._RowCount, A._ColumnCount);

            var dataC = C.Data;
            for (int i = 0; i < A._Data.Length; i++)
            {
                dataC[i] = A._Data[i] / s;
            }
            return C;
        }

        public static ComplexMatrix operator *(ComplexMatrix A, double s)
        {
            var C = new ComplexMatrix(A._RowCount, A._ColumnCount);

            var dataC = C.Data;
            for (int i = 0; i < A._Data.Length; i++)
            {
                dataC[i] = A._Data[i] * s;
            }
            return C;
        }
        public static ComplexMatrix operator *(double s, ComplexMatrix A)
        {
            return A * s;
        }


        public static explicit operator ComplexMatrix(BaseMatrix M)
        {
            var res = new ComplexMatrix(M.RowCount, M.ColumnCount,
                (i, j) => new Complex(M[i, j], 0));
            return res;
        }

        public Matrix Real()
        {
            var res = new Matrix(RowCount, ColumnCount);
            for (int i = 0; i < res.Data.Length; i++)
                res.Data[i] = Data[i].Real;
            return res;
        }
        public Matrix Imag()
        {
            var res = new Matrix(RowCount, ColumnCount);
            for (int i = 0; i < res.Data.Length; i++)
                res.Data[i] = Data[i].Imaginary;
            return res;
        }
        public Matrix Modulus()
        {
            var res = new Matrix(RowCount, ColumnCount);
            for (int i = 0; i < res.Data.Length; i++)
                res.Data[i] = Data[i].Modulus;
            return res;
        }
        public Matrix Angle()
        {
            var res = new Matrix(RowCount, ColumnCount);
            for (int i = 0; i < res.Data.Length; i++)
                res.Data[i] = Data[i].Angle;
            return res;
        }

        public ComplexMatrix Inverse()
        {
            if (this.IsSquare != true)
            {
                throw new System.ArgumentException("This is not a square matrix.");
            }
            if (ARMADILO.Enabled)
                return ARMADILO.inv(this);
            else
                throw new NotImplementedException("ComplexMatrix.Inverse: only with armadillo");
        }
        public ComplexMatrix PseudoInverse()
        {
            return SingularValueDecomposition.PseudoInverse(this);
        }
        /// <summary>
        /// Calculates the determinant of the matrix.
        /// </summary>
        /// <returns>The determinant of the matrix.</returns>
        public Complex Determinant()
        {
            if (this.IsSquare != true)
            {
                throw new System.ArgumentException("This is not a square matrix.");
            }
            if (ARMADILO.Enabled)
                return ARMADILO.det(this);
            else
                throw new NotImplementedException("ComplexMatrix.Determinant: only with armadillo");
        }
    }
}