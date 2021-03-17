using DotNumerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace DotNumerics
{
    public static class Tools
    {

        public static Vector Average<TSource>(this IEnumerable<TSource> source, Func<TSource, Vector> selector)
        {
            var v0 = selector(source.FirstOrDefault());
            var res = new Vector(v0.Type, v0.Count);
            foreach (var v in source)
            {
                var v_ = selector(v);
                for (int i = 0; i < res.Count; i++)
                    res[i] += v_[i];
            }
            res /= source.Count();
            return res;
        }
        public static Vector Average(this IEnumerable<Vector> source)
        {
            var v0 = source.FirstOrDefault();
            var res = new Vector(v0.Type, v0.Count);
            foreach (var v in source)
            {
                var v_ = v;
                for (int i = 0; i < res.Count; i++)
                    res[i] += v_[i];
            }
            res /= source.Count();
            return res;
        }
        public static Vector Max<TSource>(this IEnumerable<TSource> source, Func<TSource, Vector> selector)
        {
            var v0 = selector(source.FirstOrDefault());
            var res = new Vector(v0.Count, i => v0[i]);
            res.Type = v0.Type;
            foreach (var v in source)
            {
                var v_ = selector(v);
                for (int i = 0; i < res.Count; i++)
                    res[i] = Math.Max(res[i], v_[i]);
            }
            res /= source.Count();
            return res;
        }
        public static Vector Max(this IEnumerable<Vector> source)
        {
            var v0 = source.FirstOrDefault();
            var res = new Vector(v0.Count, i => v0[i]);
            res.Type = v0.Type;
            foreach (var v in source)
            {
                var v_ = v;
                for (int i = 0; i < res.Count; i++)
                    res[i] = Math.Max(res[i], v_[i]);
            }
            res /= source.Count();
            return res;
        }
        public static Vector Min<TSource>(this IEnumerable<TSource> source, Func<TSource, Vector> selector)
        {
            var v0 = selector(source.FirstOrDefault());
            var res = new Vector(v0.Count, i => v0[i]);
            res.Type = v0.Type;
            foreach (var v in source)
            {
                var v_ = selector(v);
                for (int i = 0; i < res.Count; i++)
                    res[i] = Math.Min(res[i], v_[i]);
            }
            res /= source.Count();
            return res;
        }
        public static Vector Min(this IEnumerable<Vector> source)
        {
            var v0 = source.FirstOrDefault();
            var res = new Vector(v0.Count, i => v0[i]);
            res.Type = v0.Type;
            foreach (var v in source)
            {
                var v_ = v;
                for (int i = 0; i < res.Count; i++)
                    res[i] = Math.Min(res[i], v_[i]);
            }
            res /= source.Count();
            return res;
        }
        public static Matrix Average<TSource>(this IEnumerable<TSource> source, Func<TSource, Matrix> selector)
        {
            var v0 = selector(source.FirstOrDefault());
            var res = new Matrix(v0.RowCount, v0.ColumnCount);
            foreach (var v in source)
            {
                var v_ = selector(v);
                for (int i = 0; i < res.Data.Length; i++)
                    res.Data[i] += v_.Data[i];
            }
            res /= source.Count();
            return res;
        }
        public static Matrix Average(this IEnumerable<Matrix> source)
        {
            var v0 = source.FirstOrDefault();
            var res = new Matrix(v0.RowCount, v0.ColumnCount);
            foreach (var v in source)
            {
                var v_ = v;
                for (int i = 0; i < res.Data.Length; i++)
                    res.Data[i] += v_.Data[i];
            }
            res /= source.Count();
            return res;
        }
        public static Matrix Max<TSource>(this IEnumerable<TSource> source, Func<TSource, Matrix> selector)
        {
            var v0 = selector(source.FirstOrDefault());
            var res = new Matrix(v0.RowCount, v0.ColumnCount, i => v0.Data[i]);
            foreach (var v in source)
            {
                var v_ = selector(v);
                for (int i = 0; i < res.Data.Length; i++)
                    res.Data[i] = Math.Max(res.Data[i], v_.Data[i]);
            }
            res /= source.Count();
            return res;
        }
        public static Matrix Max(this IEnumerable<Matrix> source)
        {
            var v0 = source.FirstOrDefault();
            var res = new Matrix(v0.RowCount, v0.ColumnCount, i => v0.Data[i]);
            foreach (var v in source)
            {
                var v_ = v;
                for (int i = 0; i < res.Data.Length; i++)
                    res.Data[i] = Math.Max(res.Data[i], v_.Data[i]);
            }
            res /= source.Count();
            return res;
        }
        public static Matrix Min<TSource>(this IEnumerable<TSource> source, Func<TSource, Matrix> selector)
        {
            var v0 = selector(source.FirstOrDefault());
            var res = new Matrix(v0.RowCount, v0.ColumnCount, i => v0.Data[i]);
            foreach (var v in source)
            {
                var v_ = selector(v);
                for (int i = 0; i < res.Data.Length; i++)
                    res.Data[i] = Math.Min(res.Data[i], v_.Data[i]);
            }
            res /= source.Count();
            return res;
        }
        public static Matrix Min(this IEnumerable<Matrix> source)
        {
            var v0 = source.FirstOrDefault();
            var res = new Matrix(v0.RowCount, v0.ColumnCount, i => v0.Data[i]);
            foreach (var v in source)
            {
                var v_ = v;
                for (int i = 0; i < res.Data.Length; i++)
                    res.Data[i] = Math.Min(res.Data[i], v_.Data[i]);
            }
            res /= source.Count();
            return res;
        }
    }
}
