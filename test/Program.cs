using System;
using System.Collections.Generic;
using System.Windows.Forms;
using DotNumerics.LinearAlgebra;
using System.Linq;

namespace test
{
    static class Program
    {
        /// <summary>
        /// The main entry point for the application.
        /// </summary>
        [STAThread]
        static void Main()
        {
            /*var K = Matrix.Import("e:\\k.txt");
            var M = Matrix.Import("e:\\m.txt");
            var C = K * 0;
            ComplexMatrix vec;
            EigenSystem.Eig(K, M, C, out vec, true, M.RowCount);*/
            //System.Numerics.Complex
            /*var res = "";
            var sw = new System.Diagnostics.Stopwatch();
            sw.Start();
            for (int h = 3; h < 145; h++)
                for (int i = 7; i < Math.Pow(2, 5); i++)
                {
                    var str = Convert.ToString(i, 2);
                    if (str.Count(c => c == '1') == 3)
                    {
                        for (int j = 0; j < str.Length; j++)
                            if (str[j] == '1')
                                res += j + "\t";
                        res += "\r\n";
                    }
                }
            Clipboard.SetText(res);
            MessageBox.Show("" + sw.ElapsedMilliseconds);
            return;*/
            Application.EnableVisualStyles();
            Application.SetCompatibleTextRenderingDefault(false);
            Application.Run(new Form1());
        }
    }
}
