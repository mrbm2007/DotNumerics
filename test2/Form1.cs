using DotNumerics;
using DotNumerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Net;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace test2
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
            checkBox1_CheckedChanged(null, null);
        }

        private void button1_Click(object sender, EventArgs e)
        {
            var e1 = new Vector(0.0, 1, 2, 3, 4, 2, 4, 2, 2, 4, 3, 5, 6, 2, 4, 6, 3).FFT();
            MessageBox.Show(e1.VectorToString());

            DotNumerics.ARMADILO.Enabled = true;
            var M = Matrix.Zeros(7, 7);
            var K = Matrix.I(7);
            var eg = EigenSystem.Eig_(K, M);
        }

        private void button_arma_Click(object sender, EventArgs e)
        {
            var check = DotNumerics.ARMADILO.ReCheck();
            var ver = -1.0;
            DateTime build = new DateTime();
            if (check)
            {
                ver = DotNumerics.ARMADILO.version();
                build = DotNumerics.ARMADILO.BuildTime();
            }
            MessageBox.Show(
                "Check:   " + check +
                "\r\nEnabled: " + DotNumerics.ARMADILO.Enabled +
                "\r\nVersion: " + ver +
                "\r\nBuild Time: " + (check ? build + "" : "")
                );
        }

        private void button3_Click(object sender, EventArgs e)
        {
            var sw1 = new System.Diagnostics.Stopwatch();
            var sw2 = new System.Diagnostics.Stopwatch();
            var n = Convert.ToInt32(textBox_n.Text);
            sw1.Start();
            var A = new ComplexMatrix(Matrix.Random(n, n, n), Matrix.Random(n, n, n * n));
            var B = Matrix.Random(n, n, 2 * n);
            //var A = new ComplexMatrix(n, n, (i, j) => new Complex(i / 7.0 + 10 * j, 9 / (2 + i)));
            //var B = new Matrix(n, 3, (i, j) => i*1.25 + 10 * j / 5.0);
            var b = B.GetColumnVector(2);
            sw1.Stop();
            //  var C = A;// * B;
            sw2.Start();
            var C = LinearEquations.Solve(A.GetReal(), A.GetImag(), b, 0.5 * b);
            //C = C.Inverse();
            textBox_res.Text = C[0].ToString1("G3") + "\r\n\r\n" + C[1].ToString1("G3") + "\r\n";// + C.Determinant();
            sw2.Stop();
            //textBox_res.Text += ddd + "   " + sw1.ElapsedMilliseconds + " , " + sw2.ElapsedMilliseconds + " , " + C[n - 1, 0] + " , " + (A.Data[n] + B.Data[n]) + "\r\n";

            if (n < 10)
            {
                var BB = new ComplexMatrix(n, n, (i, j) => (i == j ? 12 : 1) + n * i + j + 1 + new Complex(i, j));
                //var BB = new Matrix(n, n, (i, j) => n * i + j + 1);
                textBox_res.Text = BB.Determinant() + "";
                var bb = new Vector(n, i => i + 5 / (i + 1));
                var x = LinearEquations.Solve(BB, bb);
                var err = BB * x - bb;
                textBox_res.Text += "\r\n\r\n";
                textBox_res.Text += err.Norm();
                textBox_res.Text += "\r\n\r\n";
                textBox_res.Text += BB.MatrixToString("0.######");
                textBox_res.Text += "\r\n\r\n";
                var cc = BB * BB;//.Inverse();
                textBox_res.Text += cc.MatrixToString("0.######");
            }
        }

        private void checkBox1_CheckedChanged(object sender, EventArgs e)
        {
            DotNumerics.ARMADILO.Enabled = checkBox_arma.Checked;
            DotNumerics.ARMADILO.ReCheck();
            button_arma.BackColor = DotNumerics.ARMADILO.Enabled ? Color.LightGreen : Color.Orange;
        }

        private void button_clear_Click(object sender, EventArgs e)
        {
            textBox_res.Text = "";
        }

        private void button1_Click_1(object sender, EventArgs e)
        {
            var R = ARMADILO.Roots(new ComplexVector(new Complex[] {
            1,0,0,0,-25,0,0,0,0
            }));
            MessageBox.Show(R.VectorToString("G5"));
        }

        private void button2_Click(object sender, EventArgs e)
        {
            try
            {
                var A = new ComplexMatrix[]
                {
                (ComplexMatrix)new Matrix( new double[,] { { 1 } } ),
                (ComplexMatrix)new Matrix( new double[,] { { 0 } } ),
                (ComplexMatrix)new Matrix( new double[,] { { 0 } } ),
                (ComplexMatrix)new Matrix( new double[,] { { 5 } } ),
                (ComplexMatrix)new Matrix( new double[,] { { 0 } } ),
                (ComplexMatrix)new Matrix( new double[,] { { 1 } } ),
                (ComplexMatrix)new Matrix( new double[,] { { 0 } } ),
                (ComplexMatrix)new Matrix( new double[,] { { 20 } } ),
                };
                ComplexMatrix vec;
                var v = EigenSystem.Eig(A.ToList(), out vec, num: 20);
                MessageBox.Show(v.VectorToString("G5"));
            }
            catch (Exception ex)
            {
                MessageBox.Show(ex.Message);
            }
        }
    }
}
