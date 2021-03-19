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
using System.IO;
using DotNumerics.FortranLibrary;
using System.ComponentModel;
using System.Linq;

namespace DotNumerics.LinearAlgebra
{
    /// <summary>
    /// Represents a Complex Vector.
    /// </summary>
    [DebuggerDisplay(": {Type} , Length : {Length}", Name = "vector")]
    [DebuggerTypeProxy(typeof(VectorComplexDebuggerDisplay))]
    public class ComplexVector
    {

        #region Fields
        /// <summary>
        /// Los datos del vector
        /// </summary>
        [DebuggerBrowsable(DebuggerBrowsableState.Never)]
        protected Complex[] _Data;

        /// <summary>
        /// El tipo de vector.
        /// </summary>
        [DebuggerBrowsable(DebuggerBrowsableState.Never)]
        protected VectorType _Type = VectorType.Column;

        #endregion

        #region Constructor

        /// <summary>
        /// Initializes a new instance of the VectorComplex class of the given size.
        /// </summary>
        /// <param name="length">The vector length</param>
        public ComplexVector(int length) : this(VectorType.Column, length) { }

        /// <summary>
        /// Initializes a new instance of the Vector class of the given size and type.
        /// </summary>
        /// <param name="type">The vector type</param>
        /// <param name="length">length">The vector length</param>
        public ComplexVector(VectorType type, int length)
        {
            if (!BaseMatrix.AllowEmptyMatrices)
                if (length < 1) throw new System.ArgumentException("length < 1");
            this._Type = type;
            this._Data = new Complex[length];
            for (int i = 0; i < this._Data.Length; i++)
                this._Data[i] = new Complex(0, 0);
        }
        public ComplexVector(double[] real, double[] imag, VectorType type = VectorType.Column)
        {
            if (!BaseMatrix.AllowEmptyMatrices)
                if (real.Length < 1) throw new System.ArgumentException("length < 1");
            if (real.Length != imag.Length)
                throw new System.ArgumentException("ComplexVector constructor, real and imag length not equal: " + real.Length + "," + imag.Length);
            this._Type = type;
            this._Data = new Complex[real.Length];
            for (int i = 0; i < this._Data.Length; i++)
            {
                this._Data[i].Real = real[i];
                this._Data[i].Imaginary = imag[i];
            }
        }
        public ComplexVector(Vector real, Vector imag, VectorType type = VectorType.Column) : this(real.data, imag.data, type)
        {
        }

        /// <summary>
        /// Initializes a new instance of the Vector class that contains elements 
        /// copied from the specified array.
        /// </summary>
        /// <param name="data">The array whose elements are copied to the vector.</param>
        public ComplexVector(Complex[] data) : this(VectorType.Column, data) { }

        /// <summary>
        /// Initializes a new instance of the Vector class that contains elements
        /// copied from the specified array.
        /// </summary>
        /// <param name="type">The vector type</param>
        /// <param name="data">The array whose elements are copied to the vector.</param>
        public ComplexVector(VectorType type, Complex[] data)
        {
            if (!BaseMatrix.AllowEmptyMatrices)
                if (data.Length < 1) throw new System.ArgumentException("data.Length < 1");
            this._Type = type;
            this._Data = new Complex[data.Length];

            data.CopyTo(this._Data, 0);

            //for (int i = 0; i < data.Length; i++)
            //{
            //    this.MeData[i] = data[i];
            //}
        }

        public ComplexVector(int length, Func<int, Complex> filler) : this(length)
        {
            Fill(filler);
        }
        public ComplexVector(VectorType type, int length, Func<int, Complex> filler) : this(type, length)
        {
            Fill(filler);
        }
        public void Fill(Func<int, Complex> filler)
        {
            for (int i = 0; i < Count; i++)
                this[i] = filler(i);
        }

        #endregion

        #region Public Properties

        /// <summary>
        /// Los datos del vector
        /// </summary>
        [DebuggerBrowsable(DebuggerBrowsableState.Never)]
        internal Complex[] Data
        {
            get { return this._Data; }
        }


        internal double[] Data_r
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
        internal double[] Data_i
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
        /// Returns the number of elements.
        /// </summary>
        [DebuggerBrowsable(DebuggerBrowsableState.Never)]
        [EditorBrowsable(EditorBrowsableState.Never)]
        [Obsolete("Use count")]
        public int Length
        {
            get { return this._Data.Length; }
        }
        /// <summary>
        /// Returns the number of elements.
        /// </summary>
        public int Count
        {
            get { return this._Data.Length; }
        }

        public VectorType Type
        {
            get { return this._Type; }
            set { this._Type = value; }
        }

        /// <summary>
        /// Gets or sets the element at the specified index.
        /// </summary>
        /// <param name="index">The zero-based index of the element to get or set.</param>
        /// <returns>The element at the specified index.</returns>
        public virtual Complex this[int index]
        {
            get
            {
                return this._Data[index];
            }
            set
            {
                this._Data[index] = value;
            }
        }

        /// <summary>
        /// sub-vector
        /// </summary>
        /// <param name="start"></param>
        /// <param name="count">-1: to end</param>
        /// <returns></returns>
        public ComplexVector this[int start, int count]
        {
            get { return SubVector(start, count); }
            set { InsertAt(start, value, count); }
        }
        public ComplexVector this[List<int> indecis]
        {
            get { return SubVector(indecis); }
            set { InsertAt(indecis, value); }
        }
        /// <summary>
        /// sub-vector
        /// </summary> 
        public ComplexVector this[int start, int by, int to]
        {
            get { return SubVector_(start, by, to); }
        }

        public RealProp real => new RealProp(this);
        public ImagProp imag => new ImagProp(this);

        public class RealProp
        {
            ComplexVector owner;
            public RealProp(ComplexVector owner)
            {
                this.owner = owner;
            }
            public virtual double this[int i]
            {
                get { return owner._Data[i].Real; }
                set { owner._Data[i].Real = value; }
            }
        }
        public class ImagProp
        {
            ComplexVector owner;
            public ImagProp(ComplexVector owner)
            {
                this.owner = owner;
            }
            public virtual double this[int i]
            {
                get { return owner._Data[i].Imaginary; }
                set { owner._Data[i].Imaginary = value; }
            }
        }
        #endregion

        #region Operators

        /// <summary>
        /// Vector addition.
        /// </summary>
        /// <param name="A">The left side vector of the addition operator.</param>
        /// <param name="B">The right side vector of the addition operator.</param>
        /// <returns>A vector that represents the result of the addition.</returns>
        public static ComplexVector operator +(ComplexVector A, ComplexVector B)
        {
            return A.Add(B);
        }
        public static ComplexVector operator +(ComplexVector A, Vector B)
        {
            return A.Add(B);
        }
        public static ComplexVector operator +(Vector A, ComplexVector B)
        {
            return B.Add(A);
        }


        /// <summary>
        /// Unary minus.
        /// </summary>
        /// <param name="v">The vector.</param>
        /// <returns>Vector r[i] = -this[i]</returns>
        public static ComplexVector operator -(ComplexVector v)
        {
            return v.UnaryMinus();
        }


        /// <summary>
        /// Vector subtraction.
        /// </summary>
        /// <param name="A"> The left side vector of the subtraction operator.</param>
        /// <param name="B">The right side vector of the subtraction operator.</param>
        /// <returns>A vector that represents the result of the vector subtraction.</returns>
        public static ComplexVector operator -(ComplexVector A, ComplexVector B)
        {
            return A.Subtract(B);
        }
        public static ComplexVector operator -(ComplexVector A, Vector B)
        {
            return A.Subtract(B);
        }
        public static ComplexVector operator -(Vector A, ComplexVector B)
        {
            if (B.Type != A.Type || B.Count != A.Count)
                throw new ArgumentException("Vector dimensions or type are not valid: " + A.Count + " , " + B.Count);

            ComplexVector r = new ComplexVector(A.Type, A.Count);
            Complex[] rData = r.Data;
            for (int i = 0; i < rData.Length; i++)
                rData[i] = A.Data[i] - B[i];

            return r;
        }
        public static ComplexVector operator -(ComplexVector A, Complex b)
        {
            ComplexVector r = new ComplexVector(A.Type, A.Count);
            Complex[] rData = r.Data;
            for (int i = 0; i < rData.Length; i++)
                rData[i] = A.Data[i] - b;
            return r;
        }
        public static ComplexVector operator +(ComplexVector A, Complex b)
        {
            ComplexVector r = new ComplexVector(A.Type, A.Count);
            Complex[] rData = r.Data;
            for (int i = 0; i < rData.Length; i++)
                rData[i] = A.Data[i] + b;
            return r;
        }
        public static ComplexVector operator +(Complex b, ComplexVector A)
        {
            return A + b;
        }

        /// <summary>
        /// Scalar-Vector multiplication.
        /// </summary>
        /// <param name="c"> The left side complex number of the multiplication operator.</param>
        /// <param name="A">The right side vector of the multiplication operator.</param>
        /// <returns>A vector that represents the result of the multiplication.</returns>
        public static ComplexVector operator *(Complex c, ComplexVector A)
        {
            return A.Multiply(c);
        }
        public static ComplexVector operator *(ComplexMatrix A, ComplexVector b)
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
        public static ComplexVector operator *(Matrix A, ComplexVector b)
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
                    Sum += AData[i + k * ARows] * BData[k];
                CData[i] = Sum;
            }
            return C;
        }

        /// <summary>
        /// Vector-Scalar multiplication.
        /// </summary>
        /// <param name="A">The left side vector of the multiplication operator.</param>
        /// <param name="c"> The right side complex number of the multiplication operator.</param>
        /// <returns>A vector that represents the result of the multiplication.</returns>
        public static ComplexVector operator *(ComplexVector A, Complex c)
        {
            return A.Multiply(c);
        }

        ///// <summary>Vector multiplication.</summary>
        //public static Complex operator *(VectorComplex A, VectorComplex B)
        //{
        //    if (A.Type != VectorType.Row || B.Type != VectorType.Column  || B.Length != A.Length)
        //    {
        //        throw new System.ArgumentException("Vector dimensions or type are not valid.");
        //    }

        //    Complex C = new Complex(0.0, 0.0);

        //    Complex[] AData = A.Data;
        //    Complex[] BData = B.Data;

        //    for (int i = 0; i < AData.Length; i++)
        //    {
        //        C += AData[i] * BData[i];
        //    }

        //    return C;
        //}

        ///// <summary>
        ///// The dot product
        ///// </summary>
        ///// <param name="A"></param>
        ///// <returns>The dot product of A.</returns>
        //public static double DotProduct(Vector A )
        //{
        //    double C = 0.0;
        //    double[] AData = A.Data;
        //    for (int i = 0; i < AData.Length; i++)
        //    {
        //        C += AData[i] * AData[i];
        //    }
        //    return C;
        //}

        /// <summary>
        /// Transposed vector.
        /// </summary>
        /// <returns></returns>
        public ComplexVector Transpose()
        {
            ComplexVector AT = new ComplexVector(this._Data);

            if (this._Type == VectorType.Column) AT.Type = VectorType.Row;
            else AT.Type = VectorType.Column;

            return AT;
        }

        ///// <summary>Matrix- Vector multiplication.</summary>
        //public static MatrixComplex operator *(MatrixComplex A, VectorComplex B)
        //{

        //    int BRows;
        //    int BColumns;

        //    if (B.Type == VectorType.Column)
        //    {
        //        BColumns = 1;
        //        BRows = B.Length;
        //    }
        //    else
        //    {
        //        BColumns = B.Length;
        //        BRows = 1;
        //    }



        //    if (BRows != A.Columns)
        //    {
        //        throw new System.ArgumentException("Matrix dimensions are not valid.");
        //    }

        //    MatrixComplex C = new MatrixComplex(A.Rows, BColumns);

        //    Complex[] AData = A.Data;
        //    Complex[] BData = B.Data;
        //    Complex[] CData = C.Data;

        //    int ARows = A.Rows;
        //    int AColumns = A.Columns;


        //    Complex Sum;
        //    for (int j = 0; j < BColumns; j++)
        //    {
        //        for (int i = 0; i < ARows; i++)
        //        {
        //            Sum = new Complex(0.0, 0.0);
        //            for (int k = 0; k < AColumns; k++)
        //            {
        //                Sum += AData[i + k * ARows] * BData[k + j * BRows];
        //            }
        //            CData[i + j * ARows] = Sum;
        //        }
        //    }
        //    return C;
        //}


        ///// <summary>Matrix- Vector multiplication.</summary>
        //public static MatrixComplex operator *(VectorComplex A, MatrixComplex B)
        //{

        //    int ARows;
        //    int AColumns;

        //    if (A.Type == VectorType.Column)
        //    {
        //        AColumns = 1;
        //        ARows = A.Length;
        //    }
        //    else
        //    {
        //        AColumns = A.Length;
        //        ARows = 1;
        //    }



        //    if (B.Rows != AColumns)
        //    {
        //        throw new System.ArgumentException("Matrix dimensions are not valid.");
        //    }

        //    MatrixComplex C = new MatrixComplex(ARows, B.Columns);

        //    Complex[] AData = A.Data;
        //    Complex[] BData = B.Data;
        //    Complex[] CData = C.Data;

        //    int BRows = B.Rows;
        //    int BColumns = B.Columns;


        //    Complex Sum;
        //    for (int j = 0; j < BColumns; j++)
        //    {
        //        for (int i = 0; i < ARows; i++)
        //        {
        //            Sum = new Complex(0.0, 0.0);
        //            for (int k = 0; k < AColumns; k++)
        //            {
        //                Sum += AData[i + k * ARows] * BData[k + j * BRows];
        //            }
        //            CData[i + j * ARows] = Sum;
        //        }
        //    }
        //    return C;
        //}

        public static explicit operator ComplexVector(Vector M)
        {
            var res = new ComplexVector(M.Count,
                (i) => new Complex(M[i], 0));
            return res;
        }

        #endregion


        public Complex DotProduct(ComplexVector B)
        {
            var A = this;
            if (B.Count != A.Count)
                throw new System.ArgumentException("Vector dimensions must agree.");

            Complex C = 0.0;

            var AData = A.Data;
            var BData = B.Data;

            for (int i = 0; i < AData.Length; i++)
                C += AData[i] * BData[i];

            return C;
        }

        /// <summary>
        /// Implicit Vector to Matrix conversion.
        /// </summary>
        /// <param name="V">The Vector</param>
        /// <returns>The Matrix.</returns>
        public static implicit operator ComplexMatrix(ComplexVector V)
        {
            ComplexMatrix NewMatrix;
            if (V.Type == VectorType.Column)
            {
                NewMatrix = new ComplexMatrix(V.Count, 1, V.Data);
            }
            else
            {
                NewMatrix = new ComplexMatrix(1, V.Count, V.Data);
            }
            return NewMatrix;
        }

        #region Public Methods


        #region Add

        /// <summary>
        /// Add a complex number to all elements of this vector.
        /// </summary>
        /// <param name="c">The complex number.</param>
        /// <returns>
        /// VectorComplex r[i] = this[i] + c
        /// </returns>
        public ComplexVector Add(Complex c)
        {
            ComplexVector v = new ComplexVector(this._Type, this._Data.Length);
            Complex[] vData = v.Data;
            for (int i = 0; i < vData.Length; i++)
            {
                vData[i] = this._Data[i] + c;
            }

            return v;
        }

        /// <summary>
        /// In place add a scalar to all elements of this vector.
        /// </summary>
        /// <param name="c">The complex number.</param>
        public void AddInplace(Complex c)
        {
            for (int i = 0; i < this._Data.Length; i++)
            {
                this._Data[i] += c;
            }
        }


        /// <summary>
        /// Add a VectorComplex.
        /// </summary>
        /// <param name="B">The vector B.</param>
        /// <returns>
        /// VectorComplex r[i] = this[i] + B[i]
        /// </returns>
        public ComplexVector Add(ComplexVector B)
        {
            if (B.Type != this.Type || B.Count != this.Count)
                throw new ArgumentException("Vector dimensions or type are not valid.");

            ComplexVector r = new ComplexVector(this._Type, this.Count);
            Complex[] rData = r.Data;
            for (int i = 0; i < rData.Length; i++)
                rData[i] = this._Data[i] + B[i];

            return r;
        }
        public ComplexVector Add(Vector B)
        {
            if (B.Type != this.Type || B.Count != this.Count)
                throw new ArgumentException("Vector dimensions or type are not valid.");

            ComplexVector r = new ComplexVector(this._Type, this.Count);
            Complex[] rData = r.Data;
            for (int i = 0; i < rData.Length; i++)
                rData[i] = this._Data[i] + B[i];

            return r;
        }

        /// <summary>
        /// In place add a VectorComplex.
        /// </summary>
        /// <param name="B">The vector B.</param>
        public void AddInplace(ComplexVector B)
        {
            if (B.Type != this.Type || B.Count != this.Count)
            {
                throw new System.ArgumentException("Vector dimensions or type are not valid.");
            }

            for (int i = 0; i < this._Data.Length; i++)
            {
                this._Data[i] += B[i];
            }

        }

        #endregion


        #region Subtract

        /// <summary>
        /// Subtract a scalar to all elements of this vector.
        /// </summary>
        /// <param name="c">The complex number.</param>
        /// <returns>
        /// VectorComplex r[i] = this[i] - c
        /// </returns>
        public ComplexVector Subtract(Complex c)
        {
            ComplexVector v = new ComplexVector(this._Type, this.Count);
            Complex[] vData = v.Data;
            for (int i = 0; i < vData.Length; i++)
            {
                vData[i] = this._Data[i] - c;
            }

            return v;
        }

        /// <summary>
        /// In place subtract a scalar to all elements of this vector.
        /// </summary>
        /// <param name="c">The complex number.</param>
        public void SubtractInplace(Complex c)
        {
            for (int i = 0; i < this._Data.Length; i++)
            {
                this._Data[i] -= c;
            }
        }


        /// <summary>
        /// Subtract a VectorComplex.
        /// </summary>
        /// <param name="B">The vector B.</param>
        /// <returns>
        /// VectorComplex r[i] = this[i] - B[i]
        /// </returns>
        public ComplexVector Subtract(ComplexVector B)
        {
            if (B.Type != this.Type || B.Count != this.Count)
                throw new ArgumentException("Vector dimensions or type are not valid.");

            ComplexVector r = new ComplexVector(this._Type, this.Count);
            Complex[] rData = r.Data;
            for (int i = 0; i < rData.Length; i++)
                rData[i] = this._Data[i] - B[i];

            return r;
        }
        public ComplexVector Subtract(Vector B)
        {
            if (B.Type != this.Type || B.Count != this.Count)
                throw new ArgumentException("Vector dimensions or type are not valid.");

            ComplexVector r = new ComplexVector(this._Type, this.Count);
            Complex[] rData = r.Data;
            for (int i = 0; i < rData.Length; i++)
                rData[i] = this._Data[i] - B[i];

            return r;
        }

        /// <summary>
        /// In place add a VectorComplex.
        /// </summary>
        /// <param name="B">The vector B.</param>
        public void SubtractInplace(ComplexVector B)
        {
            if (B.Type != this.Type || B.Count != this.Count)
            {
                throw new System.ArgumentException("Vector dimensions or type are not valid.");
            }

            for (int i = 0; i < this._Data.Length; i++)
            {
                this._Data[i] -= B[i];
            }

        }

        #endregion


        #region Multiply


        /// <summary>
        /// Multiply a scalar to all elements of this vector.
        /// </summary>
        /// <param name="c">The complex number.</param>
        /// <returns>
        /// VectorComplex r[i] = this[i] * c
        /// </returns>
        public ComplexVector Multiply(Complex c)
        {
            ComplexVector v = new ComplexVector(this._Type, this.Count);
            Complex[] vData = v.Data;
            for (int i = 0; i < vData.Length; i++)
            {
                vData[i] = this._Data[i] * c;
            }

            return v;
        }


        /// <summary>
        /// In place multiply this vector with a scalar.
        /// </summary>
        /// <param name="scalar">The scalar </param>
        public void MultiplyInplace(Complex scalar)
        {
            for (int i = 0; i < this._Data.Length; i++)
            {
                this._Data[i] *= scalar;
            }
        }

        #endregion


        #region UnaryMinus

        /// <summary>
        /// Unary minus this vector.
        /// </summary>
        /// <returns>
        /// Vector r[i] = -this[i]
        /// </returns>
        public ComplexVector UnaryMinus()
        {
            ComplexVector v = new ComplexVector(this._Type, this.Count);
            Complex[] vData = v.Data;
            for (int i = 0; i < vData.Length; i++)
            {
                vData[i] -= this._Data[i];
            }

            return v;
        }

        /// <summary>
        /// In place unary minus of this vector.
        /// </summary>
        public void UnaryMinusInplace()
        {
            for (int i = 0; i < this._Data.Length; i++)
            {
                this._Data[i] = -this._Data[i];
            }
        }


        #endregion

        #region Conjugate

        /// <summary>
        /// Conjugate this vector.
        /// </summary>
        /// <returns>
        /// Vector r[i] = Real(this[i]) - Imaginary(this[i])
        /// </returns>
        public ComplexVector Conjugate()
        {
            ComplexVector conjVect = new ComplexVector(this._Data.Length);

            Complex[] v = conjVect.Data;
            for (int i = 0; i < v.Length; i++)
            {
                v[i] = this._Data[i].Conjugate;
            }

            return conjVect;
        }

        /// <summary>
        /// In place conjugation of this vector.
        /// </summary>
        public void ConjugateInplace()
        {
            for (int i = 0; i < this._Data.Length; i++)
            {
                this._Data[i] = this._Data[i].Conjugate;
            }
        }


        #endregion

        /// <summary>
        /// Resizes the vector
        /// Added by MRB
        /// </summary>
        /// <param name="newLength"></param> 
        public void Resize(int newLength)
        {
            Array.Resize(ref _Data, newLength);
            for (int i = Count; i < newLength - 1; i++)
                _Data[i] = new Complex(0, 0);
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
                    if (Type == VectorType.Column)
                        for (int i = 0; i < this.Count; i++)
                            A.Write(this[i].ToString(format) + "\r\n");
                    else
                        for (int i = 0; i < this.Count; i++)
                            A.Write(this[i].ToString(format) + "\t");
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
        /// Calculate the 1-norm of the vector.
        /// </summary>
        /// <returns>
        /// r = sum(abs(this[i]))
        /// </returns>
        public double Norm1()
        {
            double sum = 0;
            for (int i = 0; i < this._Data.Length; i++)
                sum += this._Data[i].Modulus;

            return sum;
        }
        public double Norm()
        {

            double norm = 0.0;

            for (int i = 0; i < this._Data.Length; i++)
                norm += Math.Pow(this._Data[i].Modulus, 2);

            norm = Math.Sqrt(norm);
            return norm;
        }

        public void InsertAt(int start, ComplexVector V, int count = -1)
        {
            if (count < 0)
                count = Math.Min(V.Count, Count - start);
            for (int i = 0; i < count; i++)
                this[i + start] = V[i];
        }
        public void InsertAt(List<int> Index, ComplexVector v1)
        {
            for (int i = 0; i < Index.Count; i++)
                this[Index[i]] = v1[i];
        }

        public ComplexVector SubVector(int start = 0, int count = -1)
        {
            if (count == -1)
                count = Count - start;
            var res = new ComplexVector(count);
            for (int i = 0; i < count; i++)
                res[i] = this[i + start];
            return res;
        }
        public ComplexVector SubVector(List<int> indecis)
        {
            var res = new ComplexVector(indecis.Count);
            for (int i = 0; i < indecis.Count; i++)
                res[i] = this[indecis[i]];
            return res;
        }
        public ComplexVector SubVector_(int from = 0, int by = 1, int to = int.MaxValue)
        {
            if (to == int.MaxValue) to = Count - 1;
            if (to < 0) to += Count - 1;
            var res_ = new List<Complex>();
            for (int i = from; i <= to; i += by)
                res_.Add(this[i]);
            return new ComplexVector(res_.ToArray());
        }
        public ComplexVector SubVector_(int from = 0, int to = int.MaxValue)
        {
            return SubVector_(from, 1, to);
        }

        ///// <summary>
        ///// Calculate the norm of the vector
        ///// </summary>
        ///// <returns>The norm</returns>
        //public double Norm()
        //{
        //    double norm = Vector.DotProduct(this);
        //    norm = Math.Sqrt(norm);
        //    return norm;
        //}


        /// <summary>
        /// Returns the equivalent string representation of the vector.
        /// </summary>
        /// <returns>The string representation of the vector.</returns>
        public string VectorToString(string format = "")
        {
            using (StringWriter writer = new StringWriter())
            {
                if (this._Type == VectorType.Column)
                {
                    for (int i = 0; i < this._Data.Length; i++)
                    {
                        writer.Write(this._Data[i].ToString(format));
                        if (i < this._Data.Length - 1) writer.WriteLine();
                    }
                }
                else if (this._Type == VectorType.Row)
                {
                    for (int i = 0; i < this._Data.Length; i++)
                    {
                        if (i < this._Data.Length - 1)
                            writer.Write(this._Data[i].ToString(format) + ", ");
                        else
                            writer.Write(this._Data[i]);
                    }
                }
                return writer.ToString();
            }
        }


        #endregion

        public double MaxAbs(int start = 0, int end = -1)
        {
            var max = 0.0;
            if (end == -1) end = Count;
            for (int i = start; i < end; i++)
                max = Math.Max(max, this[i].Modulus);
            return max;
        }
        public double MinAbs(int start = 0, int end = -1)
        {
            var min = double.MaxValue;
            if (end == -1) end = Count;
            for (int i = start; i < end; i++)
                min = Math.Min(min, this[i].Modulus);
            return min;
        }

        public int MinAbsIndex(int start = 0, int end = -1)
        {
            var min = double.MaxValue;
            var ind = 0;
            if (end == -1) end = Count;
            for (int i = start; i < end; i++)
            {
                if (min > this[i].Modulus)
                {
                    min = this[i].Modulus;
                    ind = i;
                }
            }
            return ind;
        }
        public int MaxAbsIndex(int start = 0, int end = -1)
        {
            var max = 0.0;
            var ind = 0;
            if (end == -1) end = Count;
            for (int i = start; i < end; i++)
            {
                if (max < this[i].Modulus)
                {
                    max = this[i].Modulus;
                    ind = i;
                }
            }
            return ind;
        }

        public void Swap(int i, int j)
        {
            if (i == j) return;
            var tmp = this[i];
            this[i] = this[j];
            this[j] = tmp;
        }
        /// <summary>
        /// Creates a copy of the vector.
        /// </summary>
        /// <returns>The copy of the vector.</returns>
        public ComplexVector Clone()
        {
            var NewVector = new ComplexVector(this._Type, this._Data);
            return NewVector;
        }

        /// <summary>
        /// Vector-Scalar division.
        /// </summary>
        /// <param name="A">The left side vector of the multiplication operator.</param>
        /// <param name="s"> The right side scalar of the multiplication operator.</param>
        /// <returns>A vector that represents the result of the division.</returns>
        public static ComplexVector operator /(ComplexVector A, double s)
        {
            return A.Multiply(1 / s);
        }


        internal Vector Real()
        {
            var res = new Vector(Count);
            for (int i = 0; i < Data.Length; i++)
                res.data[i] = Data[i].Real;
            return res;
        }
        internal Vector Imag()
        {
            var res = new Vector(Count);
            for (int i = 0; i < Data.Length; i++)
                res.data[i] = Data[i].Imaginary;
            return res;
        }


    }
}
