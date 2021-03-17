#region Copyright © 2009, De Santiago-Castillo JA. All rights reserved.

//Copyright © 2009 Jose Antonio De Santiago-Castillo 
//E-mail:JAntonioDeSantiago@gmail.com
//Web: www.DotNumerics.com
//
#endregion

using System;
using System.Collections.Generic;
using System.Text;
using System.Diagnostics;

//The matrix A is symmetric if it has the property A = AT, which means: 
//It has the same number of rows as it has columns; that is, it has n rows and n columns. 
//The value of every element aij on one side of the main diagonal equals its mirror 
//image aji on the other side: aij = aji for 1 <= i <= n and 1 <= j <= n. 

namespace DotNumerics.LinearAlgebra
{
    /// <summary>
    /// Represents a  symmetric matrix.
    /// </summary>
    public sealed class SymmetricMatrix : BaseMatrix
    {

        #region  Public Constructors

        /// <summary>
        /// Initializes a new instance of the SymmetricMatrix class of the given size.
        /// </summary>
        /// <param name="size">Size</param>
        public SymmetricMatrix(int size) : base(size) { }

        /// <summary>
        /// Initializes a new instance of the SymmetricMatrix class of the given size using a array
        /// </summary>
        /// <param name="size">Size</param>
        /// <param name="Data">The data</param>
        internal SymmetricMatrix(int size, double[] Data) : base(size, Data) { }

        #endregion

        #region Public Methods


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
        /// Returns the value of a element of the matrix.
        /// </summary>
        /// <param name="row">The row value (zero-based).</param>
        /// <param name="column">The column value (zero-based).</param>
        /// <returns>The matrix value at (row, column).</returns>
        public override double this[int row, int column]
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
                //aij = aji for 1 <= i <= n and 1 <= j <= n.
                this._Data[row + column * this._RowCount] =
                this._Data[column + row * this._RowCount] = value;
            }
        }



        public SymmetricMatrix Clone()
        {
            SymmetricMatrix NewMatrix = new SymmetricMatrix(this._RowCount, this._Data);
            return NewMatrix;
        }



        #region Static methods

        #region Static methods


        /// <summary>
        /// Generate a matrix with random elements
        /// </summary>
        /// <param name="size">Size</param>
        /// <returns>An m-by-n matrix with uniformly distributed
        /// random elements in <c>[0, 1)</c> interval.</returns>
        public static SymmetricMatrix Random(int size)
        {
            System.Random random = new System.Random();

            SymmetricMatrix X = new SymmetricMatrix(size);

            double[] XData = X.Data;

            for (int j = 0; j < X.ColumnCount; j++)
            {
                for (int i = 0; i < X.RowCount; i++)
                {
                    X[i, j] = random.NextDouble();
                }
            }
            return X;
        }


        /// <summary>
        /// Generate a matrix with random elements
        /// </summary>
        /// <param name="size">Size</param>
        /// <param name="Seed">
        /// A number used to calculate a starting value for the pseudo-random number
        /// sequence. If a negative number is specified, the absolute value of the number
        /// is used.
        /// </param>
        /// <returns>An m-by-n matrix with uniformly distributed
        /// random elements in <c>[0, 1)</c> interval.</returns>
        public static SymmetricMatrix Random(int size, int Seed)
        {
            System.Random random = new System.Random(Seed);

            SymmetricMatrix X = new SymmetricMatrix(size);

            double[] XData = X.Data;

            for (int j = 0; j < X.ColumnCount; j++)
            {
                for (int i = 0; i < X.RowCount; i++)
                {
                    X[i, j] = random.NextDouble();
                }
            }
            return X;
        }


        #endregion

        #endregion

        #endregion


        #region Overloading Operators

        /// <summary>
        /// Matrix addition.
        /// </summary>
        /// <param name="A">The left side matrix of the addition operator.</param>
        /// <param name="B">The right side matrix of the addition operator.</param>
        /// <returns>A matrix that represents the result of the matrix addition.</returns>
        public static SymmetricMatrix operator +(SymmetricMatrix A, SymmetricMatrix B)
        {
            if (B.RowCount != A.RowCount || B.ColumnCount != A.ColumnCount)
            {
                throw new System.ArgumentException("Matrix dimensions are not valid.");
            }

            SymmetricMatrix C = new SymmetricMatrix(A.RowCount);

            double[] AData = A.Data;
            double[] BData = B.Data;
            double[] CData = C.Data;

            for (int i = 0; i < AData.Length; i++)
            {
                CData[i] = AData[i] + BData[i];
            }

            return C;
        }

        /// <summary>
        /// Matrix subtraction.
        /// </summary>
        /// <param name="A"> The left side matrix of the subtraction operator.</param>
        /// <param name="B">The right side matrix of the subtraction operator.</param>
        /// <returns>A matrix that represents the result of the matrix subtraction.</returns>
        public static SymmetricMatrix operator -(SymmetricMatrix A, SymmetricMatrix B)
        {
            if (B.RowCount != A.RowCount || B.ColumnCount != A.ColumnCount)
            {
                throw new System.ArgumentException("Matrix dimensions are not valid.");
            }

            SymmetricMatrix C = new SymmetricMatrix(A.RowCount);

            double[] AData = A.Data;
            double[] BData = B.Data;
            double[] CData = C.Data;

            for (int i = 0; i < AData.Length; i++)
            {
                CData[i] = AData[i] - BData[i];
            }

            return C;
        }

        #region Scalar-Matrix Multiplication

        /// <summary>
        /// Scalar-Matrix multiplication.
        /// </summary>
        /// <param name="s"> The left side scalar of the multiplication operator.</param>
        /// <param name="A">The right side matrix of the multiplication operator.</param>
        /// <returns>A matrix that represents the result of the multiplication.</returns>
        public static SymmetricMatrix operator *(double s, SymmetricMatrix A)
        {
            SymmetricMatrix C = new SymmetricMatrix(A.RowCount);

            double[] AData = A.Data;
            double[] CData = C.Data;


            Matrix.MultiplicationSM(s, AData, CData);

            return C;
        }

        /// <summary>
        /// Implicit SymmetricMatrix to Matrix conversion.
        /// </summary>
        /// <param name="symmetricMatrix">The SymmetricMatrix.</param>
        /// <returns>The matrix.</returns>
        public static implicit operator Matrix(SymmetricMatrix symmetricMatrix)
        {
            Matrix NewMatrix = new Matrix(symmetricMatrix.RowCount, symmetricMatrix.ColumnCount, symmetricMatrix.Data);
            return NewMatrix;
        }

        #endregion

        #endregion

    }
}
