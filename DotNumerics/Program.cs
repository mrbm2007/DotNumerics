using System;
using System.Collections.Generic;
using System.Text;
using DotNumerics.LinearAlgebra;
using System.Diagnostics;
using System.Windows.Forms;

namespace DotNumerics
{
    static class Program
    {
        [STAThread]
        static void Main()
        {

            var frm = new Form();
            var n = 10;
            var A = Matrix.Random(n, n);
            var B = Matrix.Random(n, n);

            var sw2 = new Stopwatch();
            //BaseMatrix.thread_count = 2;
            sw2.Start();
            var C2 = -(A * B - B) + A + 2 * (7 * -B);
            sw2.Stop();

            var sw1 = new Stopwatch();
            //BaseMatrix.thread_count = 1;
            sw1.Start();
            var C1 = -(A * B - B) + A + 2 * (7 * -B);
            sw1.Stop();

            var dc = C1 - C2;
            MessageBox.Show("error: " + dc.Norm1() + "\r\nsingle: " + sw1.ElapsedMilliseconds + "\r\nmultiThread: " + sw2.ElapsedMilliseconds + "\r\nbenefit: " + (100 * sw2.ElapsedMilliseconds) / (double)sw1.ElapsedMilliseconds + "%");
            return;
            var mw = new DotNumerics.LinearAlgebra.MatrixView();
            mw.Matrix(dc);
            mw.Dock = System.Windows.Forms.DockStyle.Fill;
            frm.Controls.Add(mw);

            System.Windows.Forms.Application.Run(frm);
        }
    }
}
