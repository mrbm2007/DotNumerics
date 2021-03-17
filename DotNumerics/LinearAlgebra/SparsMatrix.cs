using System;
using System.Collections.Generic;
using System.Text;


namespace DotNumerics.LinearAlgebra
{
    /// <summary>
    /// Added By MRB
    /// </summary>
    class SparsMatrix : Matrix
    {
        public SparsMatrix(int rows, int columns)
            : base(0, 0)
        {
            
        }
        public List<int> ROW_COL;
        public List<double> VAL;
        public new double[] _Data
        {
            get { throw new Exception("''"); }
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
                var i = ROW_COL.IndexOf((int)1e6 * row + column);
                if (i >= 0)
                    return VAL[i];
                return 0;
            }
            set
            {
                if (column >= this._ColumnCount)
                {
                    throw new ArgumentException("Index was outside the bounds of the matrix (set). [" + row + "," + column + "]>[" + RowCount + "," + ColumnCount + "]");
                }
                var i = ROW_COL.IndexOf((int)1e6 * row + column);
                if (i >= 0)
                    VAL[i] = value;
                else
                {
                    VAL.Add(value);
                    ROW_COL.Add((int)1e6 * row + column);
                }
            }
        }
    }
}
