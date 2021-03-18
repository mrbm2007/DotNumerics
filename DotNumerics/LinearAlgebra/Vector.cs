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
using System.Runtime.CompilerServices;
using System.ComponentModel;
using System.Globalization;

namespace DotNumerics.LinearAlgebra
{
    /// <summary>
    /// Represents a Vector. By default a column vector
    /// </summary>
    [DebuggerDisplay(": {Type} , Length : {Length}", Name = "vector")]
    [DebuggerTypeProxy(typeof(VectorDebuggerDisplay))]
    [TypeConverter(typeof(TC_Vector))]
    public class Vector
    {

        #region Fields
        /// <summary>
        /// Los datos del vector
        /// </summary>
        [DebuggerBrowsable(DebuggerBrowsableState.Never)]
        protected double[] _Data;
        public double[] data { get { return _Data; } }
        public void SetData(double[] data)
        {
            _Data = data;
        }

        //private int MeLength = 1;

        /// <summary>
        /// El tipo de vector.
        /// </summary>
        [DebuggerBrowsable(DebuggerBrowsableState.Never)]
        protected VectorType _Type = VectorType.Column;

        #endregion

        #region Constructor
        public Vector() : this(1) { }

        /// <summary>
        /// Initializes a new instance of the Vector class of the given size.
        /// </summary>
        /// <param name="length">The vector length</param>
        public Vector(int length) : this(VectorType.Column, length) { }

        /// <summary>
        /// Initializes a new instance of the Vector class of the given size and type.
        /// </summary>
        /// <param name="type">The vector type</param>
        /// <param name="length">length">The vector length</param>
        public Vector(VectorType type, int length)
        {
            if (!BaseMatrix.AllowEmptyMatrices)
                if (length < 0) throw new System.ArgumentException("length < 0");
            this._Type = type;
            this._Data = new double[length];
        }

        /// <summary>
        /// Initializes a new instance of the Vector class that contains elements 
        /// copied from the specified array.
        /// </summary>
        /// <param name="data">The array whose elements are copied to the vector.</param>
        public Vector(params double[] data) : this(VectorType.Column, data) { }

        /// <summary>
        /// Initializes a new instance of the Vector class that contains elements
        /// copied from the specified array.
        /// </summary>
        /// <param name="type">The vector type</param>
        /// <param name="data">The array whose elements are copied to the vector.</param>
        public Vector(VectorType type, params double[] data) : this(type, true, data)
        {
        }
        /// <summary>
        /// Initializes a new instance of the Vector class that contains elements
        /// copied from the specified array.
        /// </summary>
        /// <param name="type">The vector type</param>
        /// <param name="data">The array whose elements are copied to the vector.</param>
        public Vector(VectorType type, bool MakeACopyOfData, params double[] data)
        {
            if (data.Length < 0) throw new System.ArgumentException("data.Length < 0");
            this._Type = type;

            if (MakeACopyOfData)
            {
                this._Data = new double[data.Length];
                data.CopyTo(this._Data, 0);
            }
            else
                this._Data = data;

            //for (int i = 0; i < data.Length; i++)
            //{
            //    this.MeData[i] = data[i];
            //}
        }

        public Vector(params Vector[] v)
        {
            var n = 0;
            foreach (var a in v)
                n += a.Count;
            _Data = new double[n];
            int k = 0;
            foreach (var a in v)
                foreach (var i in a._Data)
                    _Data[k++] = i;
        }


        public Vector(int length, Func<int, double> filler) : this(length)
        {
            Fill(filler);
        }
        public Vector(VectorType type, int length, Func<int, double> filler) : this(type, length)
        {
            Fill(filler);
        }
        public void Fill(Func<int, double> filler)
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
        internal double[] Data
        {
            get { return this._Data; }
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

        /// <summary>
        /// The vector type.
        /// </summary>
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
        public virtual double this[int index]
        {
            get
            {
                if (index > Count)
                    throw new System.ArgumentException("Index was outside the bounds of vector (get). [" + index + "]>[" + Count + "]");
                return this._Data[index];
            }
            set
            {
                if (index > Count)
                    throw new System.ArgumentException("Index was outside the bounds of vector (set). [" + index + "]>[" + Count + "]");
                this._Data[index] = value;
            }
        }

        /// <summary>
        /// sub-vector
        /// </summary>
        /// <param name="start"></param>
        /// <param name="count">-1: to end</param>
        /// <returns></returns>
        public Vector this[int start, int count]
        {
            get { return SubVector(start, count); }
            set { InsertAt(start, value, count); }
        }
        public Vector this[List<int> indecis]
        {
            get { return SubVector(indecis); }
            set { InsertAt(indecis, value); }
        }
        /// <summary>
        /// sub-vector
        /// </summary> 
        public Vector this[int start, int by, int to]
        {
            get { return SubVector_(start, by, to); }
        }

        public double Fisrt
        {
            get { return _Data[0]; }
            set { _Data[0] = value; }
        }
        public double Last
        {
            get { return _Data[_Data.Length - 1]; }
            set { _Data[_Data.Length - 1] = value; }
        }

        public double x
        {
            get { return this.data[0]; }
            set { this.data[0] = value; }
        }
        public double y
        {
            get { return this.data[1]; }
            set { this.data[1] = value; }
        }
        public double z
        {
            get { return this.data[2]; }
            set { this.data[2] = value; }
        }

        #endregion

        #region Operators


        /// <summary>
        /// Vector addition.
        /// </summary>
        /// <param name="A">The left side vector of the addition operator.</param>
        /// <param name="B">The right side vector of the addition operator.</param>
        /// <returns>A vector that represents the result of the addition.</returns>
        public static Vector operator +(Vector A, Vector B)
        {
            return A.Add(B);
        }


        /// <summary>
        /// Unary minus.
        /// </summary>
        /// <param name="v">The vector.</param>
        /// <returns>Vector r[i] = -this[i]</returns>
        public static Vector operator -(Vector v)
        {
            return v.UnaryMinus();
        }


        /// <summary>
        /// Vector subtraction.
        /// </summary>
        /// <param name="A"> The left side vector of the subtraction operator.</param>
        /// <param name="B">The right side vector of the subtraction operator.</param>
        /// <returns>A vector that represents the result of the vector subtraction.</returns>
        public static Vector operator -(Vector A, Vector B)
        {
            return A.Subtract(B);
        }


        /// <summary>
        /// Scalar-Vector multiplication.
        /// </summary>
        /// <param name="s"> The left side scalar of the multiplication operator.</param>
        /// <param name="A">The right side vector of the multiplication operator.</param>
        /// <returns>A vector that represents the result of the multiplication.</returns>
        public static Vector operator *(double s, Vector A)
        {
            return A.Multiply(s);
        }

        /// <summary>
        /// Vector-Scalar multiplication.
        /// </summary>
        /// <param name="A">The left side vector of the multiplication operator.</param>
        /// <param name="s"> The right side scalar of the multiplication operator.</param>
        /// <returns>A vector that represents the result of the multiplication.</returns>
        public static Vector operator *(Vector A, double s)
        {
            return A.Multiply(s);
        }

        /// <summary>
        /// Vector-Scalar division.
        /// </summary>
        /// <param name="A">The left side vector of the multiplication operator.</param>
        /// <param name="s"> The right side scalar of the multiplication operator.</param>
        /// <returns>A vector that represents the result of the division.</returns>
        public static Vector operator /(Vector A, double s)
        {
            return A.Multiply(1 / s);
        }
        /// <summary>
        /// Added By MRB
        /// </summary>
        /// <param name="A">The left side vector of the / operator.</param>
        /// <param name="s"> The right side scalar of the / operator.</param>
        /// <returns>A vector that represents the result of the A+s.</returns>
        public static Vector operator +(Vector A, double s)
        {
            var res = A.Clone();
            for (int i = 0; i < A.Count; i++)
                res.data[i] += s;
            return res;
        }
        /// <summary>
        /// Added By MRB
        /// </summary>
        /// <param name="A">The left side vector of the / operator.</param>
        /// <param name="s"> The right side scalar of the / operator.</param>
        /// <returns>A vector that represents the result of the A+s.</returns>
        public static Vector operator +(double s, Vector A)
        {
            var res = A.Clone();
            for (int i = 0; i < A.Count; i++)
                res.data[i] += s;
            return res;
        }
        /// <summary>
        /// Added By MRB
        /// </summary>
        /// <param name="A">The left side vector of the * operator.</param>
        /// <param name="B"> The right side scalar of the * operator.</param>
        /// <returns>A matrix that represents the result of the A*B.</returns>
        public static Matrix operator *(Vector A, Vector B)
        {
            return A.ToMatrix() * B.ToMatrix();
        }
        /// <summary>
        /// Added By MRB
        /// </summary>
        /// <param name="A">The left side vector of the / operator.</param>
        /// <param name="s"> The right side scalar of the / operator.</param>
        /// <returns>A vector that represents the result of the A-s.</returns>
        public static Vector operator -(Vector A, double s)
        {
            var res = A.Clone();
            for (int i = 0; i < A.Count; i++)
                res.data[i] -= s;
            return res;
        }
        /// <summary>
        /// Added By MRB
        /// </summary>
        /// <param name="A">The left side vector of the / operator.</param>
        /// <param name="s"> The right side scalar of the / operator.</param>
        /// <returns>A vector that represents the result of the s-A.</returns>
        public static Vector operator -(double s, Vector A)
        {
            var res = A.Clone();
            for (int i = 0; i < A.Count; i++)
                res.data[i] = s - res.data[i];
            return res;
        }

        /// <summary>
        /// Added by MRB
        /// </summary>
        /// <param name="i"></param>
        /// <param name="j"></param>
        public void Swap(int i, int j)
        {
            var tmp = this[i];
            this.data[i] = this[j];
            this.data[j] = tmp;
        }

        /// <summary>
        /// Added by MRB
        /// </summary>
        /// <param name="A">a row vector</param>
        /// <param name="M"></param>
        /// <returns></returns>
        public static Vector operator *(Vector A, Matrix M)
        {
            var Res = new Vector(M.ColumnCount);
            for (int c = 0; c < M.ColumnCount; c++)
            {
                double sum = 0;
                for (int r = 0; r < M.RowCount; r++)
                    sum += A.data[r] * M[r, c];
                Res.data[c] = sum;
            }
            return Res;
        }

        /// <summary>
        /// Added by MRB
        /// </summary>
        /// <param name="M"></param>
        /// <param name="A">A column vector</param>
        /// <returns></returns>
        public static Vector operator *(Matrix M, Vector A)
        {
            if (M.ColumnCount != A.Count)
                throw new Exception("M x v dimension error: " + M.RowCount + "x" + M.ColumnCount + " , " + A.Count);
            if (ARMADILO.Enabled)
                return ARMADILO.mv_multiply(M, A);

            var Res = new Vector(M.RowCount);
            for (int r = 0; r < M.RowCount; r++)
            {
                double sum = 0;
                for (int c = 0; c < M.ColumnCount; c++)
                    sum += M[r, c] * A.data[c];
                Res.data[r] = sum;
            }
            return Res;
        }

        /// <summary>
        /// Returns Zero vector with Length n
        /// Added by MRB
        /// </summary>
        /// <param name="n"></param>
        /// <returns></returns>
        public static Vector Zeros(int n)
        {
            var Res = new Vector(n);
            for (int i = 0; i < n; i++)
                Res.data[i] = 0;
            return Res;
        }

        /// <summary>
        /// Returns One vector with Length n
        /// Added by MRB
        /// </summary>
        /// <param name="n"></param>
        /// <returns></returns>
        public static Vector Ones(int n)
        {
            var Res = new Vector(n);
            for (int i = 0; i < n; i++)
                Res.data[i] = 1;
            return Res;
        }

        public void Fill(double x = 0)
        {
            for (int i = 0; i < Count; i++)
                data[i] = x;
        }


        ///// <summary>
        ///// Vector - Vector multiplication. 
        ///// Row Vector * Column Vector: Inner product.
        ///// Column Vector * Row Vector: Outer product.
        ///// </summary>
        ///// <param name="A"> The left side vector of the multiplication operator.</param>
        ///// <param name="B">The right side vector of the multiplication operator.</param>
        ///// <returns>A value that represents the result of the vector multiplication.</returns>
        ///// <remarks>
        ///// The dot product is the result of multiplying all the components of two vectors together and adding the results.
        ///// </remarks>
        //public static Matrix operator *(Vector A, Vector B)
        //{
        //    Matrix matrixA = A;
        //    Matrix matrixB = B;

        //    return  matrixA * matrixB;
        //}

        /// <summary>
        /// Dot product or scalar product.
        /// </summary>
        /// <param name="A">The left side vector of the operator.</param>
        /// <param name="B">The right side vector of the operator.</param>
        /// <remarks>
        /// The dot product is the result of multiplying all the components of two vectors together and adding the results, res= Sum(A[i]*B[i]).
        /// </remarks>
        /// <returns>The dot product = Sum(A[i]*B[i])</returns>
        /// 
#if NET45
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
        public static double DotProduct(Vector A, Vector B)
        {
            if (B.Count != A.Count)
                throw new System.ArgumentException("Vector dimensions must agree.");

            double C = 0.0;

            var AData = A.Data;
            var BData = B.Data;

            for (int i = 0; i < AData.Length; i++)
                C += AData[i] * BData[i];

            return C;
        }

        /// <summary>
        /// Dot product of this vector with another vector.
        /// </summary>
        /// <param name="B">The other vector.</param>
        /// <remarks>
        /// The dot product is the result of multiplying all the components of two vectors together and adding the results, res= Sum(A[i]*B[i]).
        /// </remarks>
        /// <returns>r = Sum(this[i]*B[i])</returns>
#if NET45
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
        public double DotProduct(Vector B)
        {
            return Vector.DotProduct(this, B);
        }
        public double Dot(Vector B) { return DotProduct(B); }

        /// <summary>
        /// Added by MRB
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <returns></returns>
#if NET45
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
        public static Vector Cross(Vector A, Vector B)
        {
            var C = new Vector(3);
            C[0] = A.data[1] * B.data[2] - A.data[2] * B.data[1];
            C[1] = A.data[2] * B.data[0] - A.data[0] * B.data[2];
            C[2] = A.data[0] * B.data[1] - A.data[1] * B.data[0];
            return C;
        }
        /// <summary>
        /// Added by MRB
        /// </summary>
        /// <param name="A"></param>
        /// <returns></returns>
#if NET45
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
        public Vector Cross(Vector A)
        {
            return Cross(this, A);
        }

        /// <summary>
        /// Transposed vector.
        /// </summary>
        /// <returns>The transposed vector.</returns>
        /// <remarks>
        /// Transposition turns a row vector into a column vector ( Or a column vector into a row vector).
        /// </remarks>
        public Vector TransposeInPlace()
        {
            if (this._Type == VectorType.Column) this.Type = VectorType.Row;
            else this.Type = VectorType.Column;
            return this;
        }

        /// <summary>
        /// Transposed vector.
        /// </summary>
        /// <returns>The transposed vector.</returns>
        /// <remarks>
        /// Transposition turns a row vector into a column vector ( Or a column vector into a row vector).
        /// </remarks>
        public Vector Transpose()
        {
            Vector AT = new Vector(this._Data);
            AT.TransposeInPlace();
            return AT;
        }
        /// <summary>
        /// Transposed vector.
        /// </summary>
        /// <returns>The transposed vector.</returns>
        /// <remarks>
        /// Transposition turns a row vector into a column vector ( Or a column vector into a row vector).
        /// </remarks>
        public Vector T()
        {
            return Transpose();
        }

        #region  Vector  And matrix Operations

        /// <summary>
        /// Matrix- Vector multiplication.
        /// </summary>
        /// <param name="M"> The left side matrix of the multiplication operator.</param>
        /// <param name="A">The right side vector of the multiplication operator.</param>
        /// <returns>A matrix that represents the result of the matrix multiplication.</returns>
        public static Matrix operator *(BaseMatrix M, Vector A)
        {
            if (M.ColumnCount != A.Count)
                throw new Exception("M x v dimension error: " + M.RowCount + "x" + M.ColumnCount + " , " + A.Count);
            if (ARMADILO.Enabled && M is Matrix)
                return ARMADILO.mv_multiply(M as Matrix, A);

            int BRows;
            int BColumns;

            if (A.Type == VectorType.Column)
            {
                BColumns = 1;
                BRows = A.Count;
            }
            else
            {
                BColumns = A.Count;
                BRows = 1;
            }


            Matrix C = new Matrix(M.RowCount, BColumns);

            double[] AData = M.Data;
            double[] BData = A.Data;
            double[] CData = C.Data;

            int ARows = M.RowCount;
            int AColumns = M.ColumnCount;


            double Sum = 0.0;
            for (int j = 0; j < BColumns; j++)
            {
                for (int i = 0; i < ARows; i++)
                {
                    Sum = 0.0;
                    for (int k = 0; k < AColumns; k++)
                    {
                        Sum += AData[i + k * ARows] * BData[k + j * BRows];
                    }
                    CData[i + j * ARows] = Sum;
                }
            }
            return C;
        }


        ///// <summary>
        ///// Vector-Matrix multiplication.
        ///// </summary>
        ///// <param name="A"> The left side vector of the multiplication operator.</param>
        ///// <param name="B">The right side matrix of the multiplication operator.</param>
        ///// <returns>A matrix that represents the result of the matrix multiplication.</returns>
        //public static Matrix operator *(Vector A, BaseMatrix B)
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

        //    Matrix C = new Matrix(ARows, B.Columns);

        //    double[] AData = A.Data;
        //    double[] BData = B.Data;
        //    double[] CData = C.Data;

        //    int BRows = B.Rows;
        //    int BColumns = B.Columns;


        //    double Sum = 0.0;
        //    for (int j = 0; j < BColumns; j++)
        //    {
        //        for (int i = 0; i < ARows; i++)
        //        {
        //            Sum = 0.0;
        //            for (int k = 0; k < AColumns; k++)
        //            {
        //                Sum += AData[i + k * ARows] * BData[k + j * BRows];
        //            }
        //            CData[i + j * ARows] = Sum;
        //        }
        //    }
        //    return C;
        //}

        #endregion

        #endregion

        /// <summary>
        /// Implicit Vector to Matrix conversion.
        /// </summary>
        /// <param name="V">The Vector</param>
        /// <returns>The Matrix.</returns>
        public static implicit operator Matrix(Vector V)
        {
            Matrix NewMatrix;
            if (V.Type == VectorType.Column)
                NewMatrix = new Matrix(V.Count, 1, V.Data);
            else
                NewMatrix = new Matrix(1, V.Count, V.Data);
            return NewMatrix;
        }
        public Matrix ToMatrix(bool MakeACopyOfData = false)
        {
            if (Type == VectorType.Column)
                return new Matrix(Count, 1, Data, MakeACopyOfData);
            else
                return new Matrix(1, Count, Data, MakeACopyOfData);
        }

        #region Public Methods
        public override string ToString()
        {
            if (Count < 10)
                return ToString1();
            return base.ToString();
        }
        public string ToString1(string format = "", bool simple = true, int digits = 4)
        {
            var res = "";
            if (!simple)
                res += "Length: " + Count + "\r\n---------\r\n";
            int max = 0;
            for (int r = 0; r < Count; r++)
            {
                if (format == "")
                {
                    var s = this[r].ToString("G" + digits);
                    max = Math.Max(max, s.Length);
                    res += '>' + s.PadLeft(digits + 6, ' ');
                }
                else
                    res += (this[r] >= 0 ? " " : "") + this[r].ToString(format);
                if (Type == VectorType.Column)
                {
                    if (!simple)
                        res += "  ;\r\n";
                    else
                        res += "\r\n";
                }
                else
                {
                    if (!simple)
                        res += ", ";
                    else
                        res += "  ";
                }
            }
            if (format == "")
            {
                var s = ">".PadRight(digits + 6 - max + 1, ' ');
                res = res.Replace(s, "");
            }
            return res;
        }

        #region To Array

        /// <summary>
        ///  Copies the elements of this vector to a new array.
        /// </summary>
        /// <returns>An array containing copies of the elements of this vector.</returns>
        public double[] ToArray()
        {

            double[] VectData = new double[this.Data.Length];

            this._Data.CopyTo(VectData, 0);

            return VectData;
        }

        #endregion


        #region Add

        /// <summary>
        /// Add a scalar to all elements of this vector.
        /// </summary>
        /// <param name="s">The scalar.</param>
        /// <returns>
        /// Vector r[i] = this[i] + s
        /// </returns>
        public Vector Add(double s)
        {
            Vector v = new Vector(this._Type, this._Data.Length);
            double[] vData = v.Data;
            for (int i = 0; i < vData.Length; i++)
            {
                vData[i] = this._Data[i] + s;
            }

            return v;
        }

        public void AddMember(double s)
        {
            Array.Resize(ref _Data, Count + 1);
            _Data[Count - 1] = s;
        }

        /// <summary>
        /// In place add a scalar to all elements of this vector.
        /// </summary>
        /// <param name="s">The scalar.</param>
        public void AddInplace(double s)
        {
            for (int i = 0; i < this._Data.Length; i++)
            {
                this._Data[i] += s;
            }
        }


        /// <summary>
        /// Add a Vector.
        /// </summary>
        /// <param name="B">The vector B.</param>
        /// <returns>
        /// Vector r[i] = this[i] + B[i]
        /// </returns>
        public Vector Add(Vector B)
        {
            if (B.Type != this.Type || B.Count != this.Count)
            {
                throw new System.ArgumentException("Vector dimensions or type are not valid.");
            }

            Vector r = new Vector(this._Type, this.Count);
            double[] rData = r.Data;
            for (int i = 0; i < rData.Length; i++)
                rData[i] = this._Data[i] + B.data[i];

            return r;
        }

        /// <summary>
        /// In place add a Vector.
        /// This[i] += B[i]
        /// </summary>
        /// <param name="B">The vector B.</param>
        public Vector AddInplace(Vector B)
        {
            if (B.Type != this.Type || B.Count != this.Count)
                throw new System.ArgumentException("Vector dimensions or type are not valid.");

            for (int i = 0; i < this._Data.Length; i++)
                this._Data[i] += B.data[i];
            return this;
        }

        #endregion


        #region Subtract

        /// <summary>
        /// Subtract a scalar to all elements of this vector.
        /// </summary>
        /// <param name="s">The scalar.</param>
        /// <returns>
        /// Vector r[i] = this[i] - s
        /// </returns>
        public Vector Subtract(double s)
        {
            Vector v = new Vector(this._Type, this.Count);
            double[] vData = v.Data;
            for (int i = 0; i < vData.Length; i++)
            {
                vData[i] = this._Data[i] - s;
            }

            return v;
        }

        /// <summary>
        /// In place subtract a scalar to all elements of this vector.
        /// </summary>
        /// <param name="s">The scalar.</param>
        public void SubtractInplace(double s)
        {
            for (int i = 0; i < this._Data.Length; i++)
            {
                this._Data[i] -= s;
            }
        }


        /// <summary>
        /// Subtract a Vector.
        /// </summary>
        /// <param name="B">The vector B.</param>
        /// <returns>
        /// Vector r[i] = this[i] - B[i]
        /// </returns>
        public Vector Subtract(Vector B)
        {
            if (B.Type != this.Type || B.Count != this.Count)
            {
                throw new System.ArgumentException("Vector dimensions or type are not valid.");
            }

            Vector r = new Vector(this._Type, this.Count);
            double[] rData = r.Data;
            for (int i = 0; i < rData.Length; i++)
            {
                rData[i] = this._Data[i] - B[i];
            }

            return r;
        }

        /// <summary>
        /// In place add a Vector.
        /// </summary>
        /// <param name="B">The vector B.</param>
        public void SubtractInplace(Vector B)
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
        /// <param name="s">The scalar.</param>
        /// <returns>
        /// Vector r[i] = this[i] * s
        /// </returns>
        public Vector Multiply(double s)
        {
            Vector v = new Vector(this._Type, this.Count);
            double[] vData = v.Data;
            for (int i = 0; i < vData.Length; i++)
                vData[i] = this._Data[i] * s;

            return v;
        }


        /// <summary>
        /// In place multiply this vector with a scalar.
        /// </summary>
        /// <param name="scalar">The scalar </param>
        public Vector MultiplyInplace(double scalar)
        {
            for (int i = 0; i < this._Data.Length; i++)
                this._Data[i] *= scalar;
            return this;
        }

        #endregion

        #region UnaryMinus

        /// <summary>
        /// Unary minus.
        /// </summary>
        /// <returns>
        /// Vector r[i] = -this[i]
        /// </returns>
        public Vector UnaryMinus()
        {
            Vector v = new Vector(this._Type, this.Count);
            double[] vData = v.Data;
            System.Threading.Tasks.Parallel.For(0, vData.Length, i =>
            {
                vData[i] -= this._Data[i];
            });
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


        #region Norm



        /// <summary>
        /// Calculate the norm of the vector (The 2-norm of the vector). 
        /// </summary>
        /// <remarks>
        /// The 2-norm of a vector is the square root of the sum of squares of the vector coefficients.
        /// res = sum(u[i]^2)
        /// </remarks>
        /// <returns>The norm</returns>
        public double Norm()
        {

            double norm = 0.0;

            for (int i = 0; i < this._Data.Length; i++)
            {
                norm += this._Data[i] * this._Data[i];
            }

            norm = Math.Sqrt(norm);
            return norm;
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
            {
                sum += Math.Abs(this._Data[i]);
            }

            return sum;
        }


        /// <summary>
        /// Calculate the p-Norm.
        /// </summary>
        /// <returns>
        /// res = Sum(abs(u[i])^p))^(1/p)
        /// </returns>
        public double NormP(int p)
        {
            if (p < 1)
            {
                throw new ArgumentOutOfRangeException("p, p < 1");
            }

            if (1 == p)
            {
                return Norm1();
            }

            if (2 == p)
            {
                return Norm();
            }

            double sum = 0;
            for (int i = 0; i < this._Data.Length; i++)
            {
                sum += Math.Pow(Math.Abs(this._Data[i]), p);
            }

            return Math.Pow(sum, 1.0 / p);
        }


        /// <summary>
        /// Infinity-Norm.
        /// </summary>
        /// <returns>
        /// res = max(abs(u[i]))
        /// </returns>
        public double NormInf()
        {
            double max = 0;
            for (int i = 0; i < this._Data.Length; i++)
            {
                max = Math.Max(max, Math.Abs(this._Data[i]));
            }

            return max;
        }


        #endregion

        /// <summary>
        /// Resizes the vector, keeps current data
        /// Added by MRB
        /// </summary> 
        /// <param name="newLength"></param>
        public void Resize(int newLength)
        {
            Array.Resize(ref _Data, newLength);
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
        public static Vector Import(string fileName)
        {
            var str = File.ReadAllText(fileName);
            var L = str.Split(new char[] { '\t', ',', ' ', '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
            var a = new List<double>();
            foreach (var c in L)
                a.Add(Convert.ToDouble(c));
            var res = new Vector(a.Count) { Type = VectorType.Row };
            var p1 = str.IndexOf('\n');
            if (p1 > 0) p1 = str.IndexOf('\n', p1 + 1);
            if (p1 > 0) p1 = str.IndexOf('\n', p1 + 1);
            if (p1 > 0) res.Type = VectorType.Column;
            for (int i = 0; i < res.Count; i++)
                res[i] = a[i];
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
            var sep = Type == VectorType.Column ? "\r\n" : "\t";
            using (var f = new StreamWriter(fileName))
            {
                for (int i = 0; i < data.Length; i++)
                    f.Write(data[i].ToString(format) + (i < data.Length - 1 ? sep : ""));
            }
        }
        /// <summary>
        /// Added by MRB, n will be normalized , Rodrigues' rotation formula
        /// </summary>
        /// <param name="n">can be un-normalized</param>
        /// <param name="theta">radian</param>
        /// <returns></returns>
        public Vector RotateAbout(Vector n, double theta)
        {
            var k = n.Normalize();
            if (k.Norm() == 0)
                k = new Vector(1e-15, 0, 0);
            var c = Math.Cos(theta);
            var s = Math.Sin(theta);
            return
                 (k.Cross(this)).MultiplyInplace(s)
                + (this * c).AddInplace(
                  k.MultiplyInplace(k.DotProduct(this) * (1 - c)));
        }
        /// <summary>
        /// Added by MRB, n will not be normalized 
        /// </summary>
        /// <param name="n">can be un-normalized</param>
        /// <param name="theta">radian</param>
        /// <returns></returns>
        public Vector RotateAbout2(Vector n, double theta = 1)
        {
            var n_ = n.Norm();
            if (n_ == 0) n_ = 1e-15;
            theta *= n_;
            var k = n / n_;
            var c = Math.Cos(theta);
            var s = Math.Sin(theta);
            return this * c
                + (k.Cross(this)).MultiplyInplace(s)
                + k.MultiplyInplace(k.DotProduct(this) * (1 - c));
        }

        public double Max()
        {
            if (Count == 0)
                return Double.NaN;
            var max = this[0];
            if (_Data.Length > 1000)
            {
                int N = 4;
                var max_ = new double[N];
                System.Threading.Tasks.Parallel.For(0, N, i =>
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
        public double AbsMax()
        {
            return ABSMax();
        }
        public double ABSMax()
        {
            if (Count == 0)
                return Double.NaN;
            var max = 0.0;
            if (_Data.Length > 1000)
            {
                int N = 4;
                var max_ = new double[N];
                System.Threading.Tasks.Parallel.For(0, N, i =>
                {
                    for (int j = i; j < _Data.Length; j += N)
                        max_[i] = Math.Max(max_[i], Math.Abs(_Data[j]));
                });
                for (int i = 0; i < N; i++)
                    max = Math.Max(max, max_[i]);
            }
            else
                foreach (var d in Data)
                    max = Math.Max(max, Math.Abs(d));
            return max;
        }
        public double Min()
        {
            if (Count == 0)
                return Double.NaN;
            var min = this[0];
            if (_Data.Length > 1000)
            {
                int N = 4;
                var min_ = new double[N];
                System.Threading.Tasks.Parallel.For(0, N, i =>
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

        public Vector Abs()
        {
            var res = new Vector(Count);
            for (int i = 0; i < Count; i++)
                res[i] = Math.Abs(this[i]);
            return res;
        }
        public Vector Sqrt()
        {
            var res = new Vector(Count);
            for (int i = 0; i < Count; i++)
                res[i] = Math.Sqrt(this[i]);
            return res;
        }
        public Vector Sign()
        {
            var res = new Vector(Count);
            for (int i = 0; i < Count; i++)
                res[i] = Math.Sign(this[i]);
            return res;
        }
        public double Sum()
        {
            var res = 0.0;
            for (int i = 0; i < Count; i++)
                res += this[i];
            return res;
        }

        public string ToString_3col(string format = "")
        {
            var res = "";
            for (int i = 0; i < data.Length - 2; i += 3)
                res += Data[i].ToString(format).PadRight(15) + " " +
                    Data[i + 1].ToString(format).PadRight(15) + " " +
                    Data[i + 2].ToString(format).PadRight(15) + "\r\n";
            return res;
        }
        public string ToString_6col(string format = "")
        {
            var res = "";
            for (int i = 0; i < data.Length - 5; i += 6)
                res += Data[i].ToString(format).PadRight(15) + " " +
                    Data[i + 1].ToString(format).PadRight(15) + " " +
                    Data[i + 2].ToString(format).PadRight(15) + " " +
                    Data[i + 3].ToString(format).PadRight(15) + " " +
                    Data[i + 4].ToString(format).PadRight(15) + " " +
                    Data[i + 5].ToString(format).PadRight(15) + "\r\n";
            return res;
        }
        public Vector Round(int digits = 0)
        {
            var res = new Vector(Count);
            for (int i = 0; i < Data.Length; i++)
                res.Data[i] = Math.Round(Data[i], digits);
            return res;
        }

        /// <summary>
        /// Normalizes this vector to a unit vector with respect to the Eucliden 2-Norm.
        /// </summary>
#if NET45
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
        public Vector Normalize()
        {
            double norm = Norm();
            Vector normalized = this.Clone();
            if (norm < 1E-13)
            {
                return normalized;
            }
            normalized.MultiplyInplace(1.0 / norm);
            return normalized;
        }
        /// <summary>
        /// Added by MRB
        /// </summary>
#if NET45
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
        public Vector NormalizeInPlace()
        {
            double norm = Norm();
            if (norm > 0)
                MultiplyInplace(1.0 / norm);
            return this;
        }

        /// <summary>
        /// Creates a copy of the vector.
        /// </summary>
        /// <returns>The copy of the vector.</returns>
        public Vector Clone()
        {
            Vector NewVector = new Vector(this._Type, this._Data);
            return NewVector;
        }



        public static explicit operator ComplexMatrix(Vector v)
        {
            var v2 = new ComplexVector(v.Count);
            for (int i = 0; i < v.Count; i++)
                v2.Data[i].Real = v[i];
            return v2;
        }

        /// <summary>
        /// Returns the equivalent string representation of the vector.
        /// </summary>
        /// <returns>The string representation of the vector.</returns>
        public string VectorToString()
        {
            using (StringWriter writer = new StringWriter())
            {
                if (this._Type == VectorType.Column)
                {
                    for (int i = 0; i < this._Data.Length; i++)
                    {
                        writer.Write(this._Data[i]);
                        if (i < this._Data.Length - 1) writer.WriteLine();
                    }
                }
                else if (this._Type == VectorType.Row)
                {
                    for (int i = 0; i < this._Data.Length; i++)
                    {
                        if (i < this._Data.Length - 1)
                            writer.Write(this._Data[i] + ", ");
                        else
                            writer.Write(this._Data[i]);
                    }
                }
                return writer.ToString();
            }
        }

        public void InsertAt(int start, Vector V, int count = -1)
        {
            if (count < 0)
                count = Math.Min(V.Count, Count - start);
            for (int i = 0; i < count; i++)
                this[i + start] = V[i];
        }
        public void InsertAt(List<int> Index, Vector v1)
        {
            for (int i = 0; i < Index.Count; i++)
                this[Index[i]] = v1[i];
        }

        public Vector SubVector(int start = 0, int count = -1)
        {
            if (count == -1)
                count = Count - start;
            var res = new Vector(count);
            for (int i = 0; i < count; i++)
                res[i] = this[i + start];
            return res;
        }
        public Vector SubVector(List<int> indecis)
        {
            var res = new Vector(indecis.Count);
            for (int i = 0; i < indecis.Count; i++)
                res[i] = this[indecis[i]];
            return res;
        }
        public Vector SubVector_(int from = 0, int by = 1, int to = int.MaxValue)
        {
            if (to == int.MaxValue) to = Count - 1;
            if (to < 0) to += Count - 1;
            var res_ = new List<double>();
            for (int i = from; i <= to; i += by)
                res_.Add(this[i]);
            return new Vector(res_.ToArray());
        }
        public Vector SubVector_(int from = 0, int to = int.MaxValue)
        {
            return SubVector_(from, 1, to);
        }

        /// <summary>
        /// Added by MRB
        /// returns small angle between this and v, (0..PI)
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        public double AngleWith(Vector v, Vector n = null)
        {
            double r = this.DotProduct(v) / (this.Norm() * v.Norm());
            if (r > 1)
                r = 1;
            else if (r < -1)
                r = -1;
            if (n != null)
                if (n.AngleWith(this.Cross(v)) > Math.PI / 2)
                    r = 2 * Math.PI - r;
            return Math.Acos(r);
        }


        /// <summary>
        /// Added by MRB
        /// </summary>
        /// <param name="x1"></param>
        /// <param name="x2"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public static Vector Linspace_(double x1, double x2, int n = 15, int i1 = 0)
        {
            var res = new List<double>();
            for (int i = i1; i < n; i++)
                res.Add(x1 + i * (x2 - x1) / (n - 1));
            return new Vector(res.ToArray());
        }
        /// <summary>
        /// Added by MRB
        /// </summary>
        /// <param name="x1"></param>
        /// <param name="x2"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public Vector Linspace(double x1, double x2, int n = 15, int i1 = 0)
        {
            var res = new List<double>();
            for (int i = i1; i < n; i++)
                res.Add(x1 + i * (x2 - x1) / (n - 1));
            _Data = res.ToArray();
            return this;
        }

        public static Vector Concatenate(params Vector[] V)
        {
            var res = new List<double>();
            foreach (var v in V)
                res.AddRange(v.data);
            return new Vector(res.ToArray());
        }

        public static Vector Combine(Vector v1, Vector v2)
        {
            var v = new Vector(v1.Count + v2.Count);
            v.InsertAt(0, v1);
            v.InsertAt(v1.Count, v2);
            return v;
        }
        public void Split(int FirstCount, ref Vector v1, ref Vector v2)
        {
            v1 = this[0, FirstCount];
            v2 = this[FirstCount, -1];
        }
        public void Split(List<int> Index1, ref Vector v1, ref Vector v2)
        {
            var Index2 = new List<int>();
            for (int i = 0; i < Count; i++)
                if (!Index1.Contains(i))
                    Index2.Add(i);
            v1 = v1 ?? new Vector(Index1.Count);
            v2 = v2 ?? new Vector(Index2.Count);
            for (int i = 0; i < Index1.Count; i++)
                v1[i] = this[Index1[i]];

            for (int i = 0; i < Index2.Count; i++)
                v2[i] = this[Index2[i]];
        }
        public static Vector Cobine(List<int> Index1, Vector v1, Vector v2)
        {
            var v = new Vector(v1.Count + v2.Count);
            var Index2 = new List<int>();
            for (int i = 0; i < v.Count; i++)
                if (!Index1.Contains(i))
                    Index2.Add(i);
            for (int i = 0; i < Index1.Count; i++)
                v[Index1[i]] = v1[i];

            for (int i = 0; i < Index2.Count; i++)
                v[Index2[i]] = v2[i];

            return v;
        }

        public ComplexVector FFT(int N = 0)
        {
            if (ARMADILO.Enabled)
                return ARMADILO.FFT(this, N);
            else
                throw new NotImplementedException("FFT with ARMADILO only!");
        }
        public void FFT(out Vector result_real, out Vector result_imag, int N = 0)
        {
            if (ARMADILO.Enabled)
                ARMADILO.FFT(this, out result_real, out result_imag, N);
            else
                throw new NotImplementedException("FFT with ARMADILO only!");
        }

        #endregion
    }

    public class TC_Vector : TypeConverter
    {
        public static string format = "0.######";
        public override bool CanConvertFrom(ITypeDescriptorContext context, Type sourceType)
        {
            return sourceType == typeof(string);
        }
        public override bool CanConvertTo(ITypeDescriptorContext context, Type destinationType)
        {
            return destinationType == typeof(string);
        }
        public override object ConvertTo(ITypeDescriptorContext context, CultureInfo culture, object value, Type destinationType)
        {
            if (value == null) return "?";
            var v_ = (value as Vector);
            var res = "";
            foreach (var d in v_.data)
                res += d.ToString(format) + ", ";
            res = "<" + res.Substring(0, res.Length - 2) + ">";
            return res;
        }
        public override object ConvertFrom(ITypeDescriptorContext context, CultureInfo culture, object value)
        {
            var str = value.ToString();
            var A = str.Split(new char[] { ',', ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
            var res = new Vector(A.Length);
            int i = 0;
            foreach (var a in A)
                res[i++] = Convert.ToDouble(a);
            return res;
        }
    }
}
