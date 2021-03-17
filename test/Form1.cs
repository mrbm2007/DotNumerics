using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Text;
using System.Windows.Forms;
using DotNumerics;
using DotNumerics.LinearAlgebra;

namespace test
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
            //   ARMADILO.dir = @"D:\Documents\Visual Studio 2012\Projects\DotNumeric\armadillo4cs\dll\";
            if (!ARMADILO.Check())
            {
                checkBox_ARMADILO.Checked = checkBox_ARMADILO.Enabled = false;
                Text = "Armadillo check failed";
            }
            else
            {
                ARMADILO.Enabled = true;
                Text = "Armadillo build: " + ARMADILO.BuildTime() + "";
            }
            checkBox_ARMADILO.Checked = ARMADILO.Enabled;
        }

        private void button1_Click(object sender, EventArgs e)
        {
            if (!ARMADILO.Check())
            {
                ARMADILO.ReCheck();
                if (!ARMADILO.Check())
                {
                    checkBox_ARMADILO.Checked = checkBox_ARMADILO.Enabled = false;
                    Text = "Armadillo check failed";
                }
                else
                {
                    ARMADILO.Enabled = true;
                    Text = "Armadillo build: " + ARMADILO.BuildTime() + "";
                }
            }
            ARMADILO.Enabled = checkBox_ARMADILO.Checked;

            int n = Convert.ToInt32(textBox3.Text);
            int m = Convert.ToInt32(textBox3.Text);

            var A = Matrix.Random(n, m, DateTime.Now.Millisecond);
            var B = Matrix.Random(n, m, DateTime.Now.Millisecond);

            var sw1 = new System.Diagnostics.Stopwatch();
            sw1.Start();
            var x3 = A * B;
            sw1.Stop();
            if (n + m < 20)
            {
                textBox2.Text = x3.ToString1();
            }
            textBox1.Text += sw1.ElapsedMilliseconds + " - ";
            Refresh();
            Application.DoEvents();
            ARMADILO.Enabled = false;

            //MessageBox.Show(this, "");

            var A2 = ILNumerics.ILMath.array<double>(A.Data, n, m);
            var B2 = ILNumerics.ILMath.array<double>(B.Data, n, m);
            var sw2 = new System.Diagnostics.Stopwatch();
            sw2.Start();
            var x4_ = ILNumerics.ILMath.multiply(A2, B2);
            var x4 = new Matrix(n, m);
            var d = x4_.GetArrayForRead();
            sw2.Stop();
            Array.Resize(ref d, x4.Data.Length);
            x4.SetData(d);

            if (n + m < 20)
            {
                textBox2.Text += "\r\n\r\n\r\n";
                textBox2.Text += x4_.ToString();
                textBox2.Text += "\r\n";
                textBox2.Text += "\r\n";
            }
            var dc = x4 - x3;
            dc.ElemntsDiv((x4 + x3) / 2.0);
            textBox2.Text += "Error:" + dc.MaxABS() * 100 + " %\r\n";

            textBox1.Text += sw2.ElapsedMilliseconds + " - ";
            if (sw1.ElapsedMilliseconds > 0)
                textBox1.Text += (100 * sw2.ElapsedMilliseconds / sw1.ElapsedMilliseconds) + " %";
            textBox1.Text += "\r\n";

            Refresh();
            Application.DoEvents();
            ARMADILO.Enabled = false;

            if (n + m < 20)
            {
                textBox1.Text += "\r\n\r\n\r\n";
                textBox1.Text += A.ToString1();
            }

            if (checkBox1.Checked)
            {
                A.Export4Matlab("test.m", "A");
                B.Export4Matlab("test.m", "B", true);
                x3.Export4Matlab("test.m", "x3", true);
                x4.Export4Matlab("test.m", "x4", true);
            }
        }

        private void button2_Click(object sender, EventArgs e)
        {
            MessageBox.Show("ARMADILO : " + ARMADILO.ReCheck());
        }

        private void button3_Click(object sender, EventArgs e)
        {

            var M = Matrix.Zeros(7, 7);
            EigenSystem.Eig(M, M);
        }
    }
}
